"""
memspectrum
===========

A package for Maximum Entropy Spectral Analysis.
Implements Burg's method to compute the power spectral density of a given 
time-series. It can be used with two different methods, the standard, which
provide higher stability, and a faster version. Having solved the method, 
it is possible to use the algorithm to forecast on the time-series. 

"""

import sys
import numpy as np
import warnings

#Stefano: should the sampling rate be given at initialization? Or is it fine like it is today?

#############DEBUG LINE PROFILING
try:
    from line_profiler import LineProfiler

    def do_profile(follow=[]):
        def inner(func):
            def profiled_func(*args, **kwargs):
                try:
                    profiler = LineProfiler()
                    profiler.add_function(func)
                    for f in follow:
                        profiler.add_function(f)
                    profiler.enable_by_count()
                    return func(*args, **kwargs)
                finally:
                    profiler.print_stats()
            return profiled_func
        return inner
except:
    pass

#pip install line_profiler
#add decorator @do_profile(follow=[]) before any function you need to track

#################

class optimizer:

    def __init__(self, method):
        """
        Implements various method to choose the best recursive order for Burg's Algorithm
        Avilable methods are "FPE", "OBD", "CAT", "AIC". The most representative order is 
        chosen to be the one that minimized the related function. 
        Parameters
        ----------
        method : 'str'
            Selects the method to be used to estimate the best order between "FPE",
            "OBD", "CAT", "AIC"


        """
        self.method = method

    def __call__(self, *args): #Stefano: what are args? We should specify them and call them by name...
        if self.method == 'FPE':
            return self._FPE(args[0], args[2], args[3])
        elif self.method == 'CAT':
            return self._CAT(args[0], args[2], args[3])
        elif self.method == 'OBD':
            return self._OBD(*args)
        elif self.method =='AIC':
            return self._AIC(args[0], args[2], args[3])
        elif self.method == 'Fixed':
            return self._Fixed(args[3])
        else:
            raise ValueError("{} is not a an available method! Valid choices are 'FPE', 'AIC', 'CAT', 'OBD' and 'Fixed'.".format(self.method))
    
    def _FPE(self, P, N, m):
        """
        Implements Akaike Final prediction Error to estimate the recursive 
        order 
        
        Parameters
        ----------
        P : 'np.float'
            The estimate of the variance for the white noise component.
        N : 'np.int'
            The length of the dataset.
        m : 'np.int'
            The recursive order.

        Returns
        -------
        'np.float'
            The value of FPE optimizer.

        """
        return P[-1] * (N + m + 1) / (N - m - 1)
    
    def _AIC(self, P, N, m):
        """
        Implements Akaike information criterion to estimate the recursive 
        order 
        
        Parameters
        ----------
        P : 'np.float'
            The estimate of the variance for the white noise component.
        N : 'np.int'
            The length of the dataset.
        m : 'np.int'
            The recursive order.

        Returns
        -------
        'np.float'
            The value of AIC optimizer.

        """
        return np.log(P[-1]) + (2 * m) / N
    
    def _CAT(self, P, N, m):
        """
        Implements Parzen's criterion on autoregressive transfer function
        to estimate the recursive order 
        
        Parameters
        ----------
        P : 'np.float'
            The estimate of the variance for the white noise component.
        N : 'np.int'
            The length of the dataset.
        m : 'np.int'
            The recursive order.

        Returns
        -------
        'np.float'
            The value of CAT optimizer.

        """
        if m == 0:
            return np.inf
        P = np.array(P[1:])
        k = np.linspace(1, m, m)
        PW_k = N / (N - k) * P
        return 1 / (N * PW_k.sum())- (1 / PW_k[-1])
    
    def _OBD(self, P, a_k, N, m):
        """
        Implement Rao's Optimum Bayes Decision rule to estimate the recursive
        order

        Parameters
        ----------
         P : 'np.float'
            The estimate of the variance for the white noise component.
        a_k : 'np.array'
            The values for the final prediction error coefficients for the 
            considered recursive order.
        N : 'np.int'
            The length of the dataset.
        m : 'np.int'
            The recursive order.

        Returns
        -------
        TYPE
            The value of OBD optimizer.

        """
        P_m = P[-1]
        P = np.array(P[:-1])
        return (N - m - 2)*np.log(P_m) + m*np.log(N) + np.log(P).sum() + (a_k**2).sum()
    
    def _Fixed(self, m):
        return 1./m

class MESA(object):
    """
    Class the implement the reproduces the Maximum Entropy Spectrum of a given 
    time-series. 
    
    init: data: `np.ndarray` shape (N,)
    solve(): 
    """
    def __init__(self, data):
        """ 
        Class that implements Burg method to estimate power spectral densitiy of
        time series. 
        
        Parameters
        ----------
        data: 'np.ndarray'       
            Time series with power spectral density to be computed 
        """ 
        self.data = data
        self.N    = len(self.data)
        self.P = None
        self.a_k = None #If a_k and P are None, the model is not already fitted
        self.optimization = None

    def save(self,filename):
        """
        Save the class to file (if the spectral density analysis is already performed).
        The output file can be used to load the class with method load()
        File is a 1D array with the format: [P, a_k, optimization, data]. The header holds the shapes for each array
        
        Parameters
        ----------
        filename: `string`      
            Name of the file to save the data at
        """
        if self.P is None or self.a_k is None:
            raise RuntimeError("PSD analysis is not performed yet: unable to save the model. You should call solve() before saving to file") 
        
        to_save = np.concatenate([[self.P], self.a_k, self.optimization, self.data])
        header = "(1,{},{},{})".format(len(self.a_k), len(self.optimization), len(self.data))
        print(header)
        np.savetxt(filename, to_save, header = header)
        return
        
    def load(self,filename):
        """
        Load the class from a given file. The file shall be the same format produced by save().
        
        Parameters
        ----------
        filename: `string`      
            Name of the file to load the data from
        """
        data = np.loadtxt(filename)
        if data.ndim != 1:
            raise ValueError("Wrong format for the file: unable to load the model")
        
            #reading the first line
        with open(filename) as f:
            first_line = f.readline()

        shapes = eval(first_line.translate({ord('#'): None, ord(' '): None}))
            #checking for the header
        if not isinstance(shapes, tuple):
            if len(shapes) != 4 or not np.all([isinstance(s, int) for s in shapes]):
                raise ValueError("Wrong format for the header: unable to load the model")

            #assigning values
        self.P, self.a_k, self.optimization, self.data = np.split(data, np.cumsum(shapes)[:-1])

        self.N = shapes[3]
        
        return
        
    def save_spectrum(self, filename, dt, frequencies = None):
        """
        Saves the power spectral density computed by the model to a txt file. The PSD is evaluated on a user given grid of frequencies. If None, a standard grid is used (as computed by np.fft.fftfreq).
        The spectrum is saved as a 2D array: [f, PSD(f)]
        
        Parameters
        ----------
        filename: `string`      
            Name of the file to save the PSD at

        dt: `np.float`      
            Sampling rate for the time series

        frequencies: `np.ndarray`      
            Frequencies to compute the PSD at. If None, a standard grid will be computed
                
        """
        if isinstance(frequencies, np.ndarray):
            PSD = self.spectrum(dt,frequencies)
        elif frequencies is None:
            PSD, frequencies = self.spectrum(dt)
        else:
            raise ValueError("Wrong type for input frequencies: expected \`np.ndarray\` or None, got {} instead".format(type(frequencies)))
         
        to_save  = np.column_stack([frequencies, PSD])
        np.savetxt(filename, to_save, header = "dt = {} s".format(dt))
        return
        
        
    def _spectrum(self, dt, N):
        """
        Method that compted the spectrum of the time series on sampling 
        frequencies 

        Parameters
        ----------
        dt: `np.float`      
            Sampling rate for the time series
            
        N: `np.float`       
            Length of the frequency grid

        Returns
        -------
        spectrum: `np.ndarray`   
            PSD of the model, including both positive and negative frequencies (shape (N,))
            
        freq: `np.ndarray`       
            Frequencies at which spectrum is evaluated (as provided by np.fft.fftfreq) (shape (N,))
        """
       
        den = np.fft.fft(self.a_k, n=N)
        spectrum = dt * self.P / (np.abs(den) ** 2)
        return spectrum, np.fft.fftfreq(N,dt)



    def spectrum(self, dt, frequencies = None): 
        """
        Computes the power spectral density of the model. Default returns power 
        spectral density and frequency array automatically computed by sampling theory. 
        It can also be computed on a user-given frequency grid passing proper frequency
        array
        
        Parameters: 
        ----------
        dt: 'np.float'                   
            Sampling rate for the time series 
            
        frequencies: 'np.ndarray'        
            (positive) frequencies to evaluate the spectrum at (shape (N,))

        

        Returns: 
        ----------
        if no frequency array is given 
            spectrum: 'np.ndarray'           
                Two sided power spectral density (shape = (N,))
            frequencies: 'np.ndaarray'      
                Frequencies at which power spectral density is evaluated (shape = (N,))
            
        if frequency array is given:
            spectrum: 'np.ndarray'           
                Power spectral density interpolated on desidered frequencies (shape = (N,))
    
            
        Raises: 
        ----------
            ValueError if frequencies greater then Nyquist frequencies are given 
            
        """
        f_ny = .5 / dt 
        
        if isinstance(frequencies, np.ndarray): 
            df = np.min(np.abs(np.diff(frequencies))) * 0.9
            if np.max(frequencies) > f_ny *1.01: 
                warnings.warn("Some of the required frequencies are higher than the Nyquist frequency ({} Hz): a zero PSD is returned for f>Nyquist".format(f_ny))
        
        if frequencies is None: N = self.N
        elif isinstance(frequencies, np.ndarray): N = int(2. * f_ny / df)
        else: raise ValueError("Type of frequencies not understood: expected to be None or np.ndarray but given {} insted".format(type(frequencies)))

        spec, f_spec = self._spectrum(dt, N)
        
        if frequencies is None: 
            return spec, f_spec   

        f_interp = np.interp(frequencies, f_spec[:int(N/2+0.5)], spec.real[:int(N/2+0.5)], left = 0., right = 0.)
        
        return f_interp 

        
    def solve(self,
              m = None,
              optimisation_method = "FPE",
              method              = "Fast",
              regularisation      = 1.0,
              early_stop = True ):
        """
        Computes the power spectral density of the attribute data for the class
        using standard Burg method recursive and a Faster (but less stable) version. 
        Default is Fast. 

        Parameters
        ----------
        m : 'np.int'                   
            Maximum number of recursions for the computation of the  power spectral density. 
            Default is None, that means m = 2N / log(2N)
                                 
        optimisation_method: 'str'     
            Method used to select the best recursive order. The order is chosen
            minimizing the corresponding method. 
            Available methods are "FPE", "OBD", "CAT", "AIC".
            Deafult is "FPE".   
        
        method: 'str'                  
            Can be "standard" or "Fast". Selects the algorithm  used to compute 
            the power spectral density with. Default is "Fast"
                                       
        regularisation: 'np.float'     
            Implement Tikhonov regularisation. Should be a number slightly larger than one. 
            Default is 1, which means no regularisation 
                                       
        early_stop: 'boolean'          
            Default is True. Breaks the iteration if there is no  no new global 
            maximum after 200 iterations. 
            Recommended for every optimisation method but CAT.
        

        Returns
        -------
        P: 'np.float'                  
            Variance of white noise for the associated autoregressive process 
                                       
        a_k: 'np.ndarray'             
            The coefficient used to compute the power spectral density (Shape (N,)) 
            
        optimization: 'np.ndarray'    
            The values of the chosen optimisation_method at every iteration 
            (Shape (N,))   

        """
        
        self.regularisation = regularisation
        self.early_stop = early_stop
        if m is None:
            self.mmax = int(2*self.N/np.log(2.*self.N))
        else:
            self.mmax = m
           
        if method == "Fast":
            self._method = self._FastBurg
        elif method == "Standard":
            self._method = self._Burg
        else:
            print("Method {0} unknown! Valid choices are 'Fast' and 'Standard'".format(method))
            exit(-1)
        
        self._optimizer = optimizer(optimisation_method)
        self.P, self.a_k, self.optimization = self._method()
        return self.P, self.a_k, self.optimization

    #@do_profile(follow=[])
    def _FastBurg(self):
        """
        Uses the Fast version of Burg Algorithm to compute the power spectral
        density. The order is selected by the minimization of the chosen method

        Returns
        -------
        P: 'np.float'                  
            Variance of white noise for the associated autoregressive process
                                       
        a_k: 'np.ndarray'             
            The coefficient used to compute the power spectral density (Shape (N,)) 
            
        optimization: 'np.ndarray'    
            The values of the chosen optimisation_method at every iteration 
            (Shape (N,))   
        """
        
        #Define autocorrelation
        c = np.zeros(self.mmax + 2, dtype = self.data.dtype) #here c has the same type of the data
        #FIXME: use numpy functions (Stefano: not really simple to do this.. Now it is the bottleneck of the function)
        for j in range(self.mmax + 1):
            c[j] = self.data[: self.N - j] @ self.data[j : ]
        c[0] *= self.regularisation
        #Initialize variables
        a = [np.array([1])]
        P = [c[0] / self.N]
        r = 2 * c[1]
        g = np.array([2 * c[0] - self.data[0] * self.data[0].conj() - self.data[-1] * self.data[-1].conj(), r])
        #Initialize lists to save arrays
        optimization = []
        idx = None
        old_idx = 0
        #Loop for the FastBurg Algorithm
        for i in range(self.mmax):
            #Update prediction error filter and reflection error coefficient
            k, new_a = self._updateCoefficients(a[-1], g)
            #Update variables. Check paper for indeces at j-th loop.
            r = self._updateR(i, r, c[i + 2])
            #Construct the array in two different, equivalent ways.
            DrA = self._constructDr2(i, new_a)
            #Update last coefficient
            g = self._updateG(g, k, r, new_a, DrA)
            #Append values to respective lists
            a.append(new_a)
            P.append(P[-1] * (1 - k * k.conj()))
            #Compute optimizer value for chosen method
            optimization.append( self._optimizer(P, a[-1], self.N, i + 1) )
            
                #checking if there is a minimum (every some iterations) if early_stop option is on
            if ((i % 200 == 0 and i !=0) or (i >= self.mmax-1)) and self.early_stop:
                idx = np.argmin(optimization) + 1
                if old_idx < idx:
                    old_idx = idx
                else:
                    old_idx = idx
                    break
        if not self.early_stop:
            idx = np.argmin(optimization) + 1

        return P[idx], a[idx], np.array(optimization)
   
    def _updateCoefficients(self, a, g):
        """
        Updates the forward prediction error coefficients (needed to compute 
        to compute the spectrum) from order i to order i + 1 
        i + 1

        Parameters
        ----------
        a : 'np.ndarray'
            The i'th order forward prediction error coefficients. (Shape (i,))
        g : 'np.ndarray'
            Array used to update the forward prediction error (Shape (i+1, )).

        Returns
        -------
        k : 'np.float'
            The reflection coefficient.
        aUpd : 'np.ndarray'
            The (i+1)'th order forward predicition error coefficients (Shape (i+1,)).

        """
        a = np.concatenate((a, np.zeros(1)))
        k = - (np.dot(a.conj(), g[::-1])) / (np.dot(a, g))
        aUpd = a + k * a[::-1].conj()
        return k, aUpd

    def _updateR(self, i, r, aCorr):
        rUp = np.array([2 * aCorr])
        rDown = r - self.data[: i + 1] * self.data[i + 1].conj() - \
            self.data[self.N - i - 1 : ].conj()[::-1] * self.data[self.N - i - 2]
        return np.concatenate((rUp, rDown))
        
    def _constructDr(self, i, a):
        #print(i, 'Outer')
        data1 = self.data[ : i + 2][::-1]
        data2 = self.data[self.N - i - 2 :]
        d1 = -np.outer(data1, data1.conj())
        d2 = -np.outer(data2.conj(), data2)
        return d1 + d2

    def _constructDr2(self, i, a):
        data1 = self.data[ : i + 2][::-1]
        data2 = self.data[self.N - i - 2 :]
        d1 = - data1 * np.dot(data1, a).conj()
        d2 = - data2.conj() * np.dot(data2, a.conj())
        return d1 + d2

    def _updateG(self, g, k, r, a, Dra):
        gUp = g + (k * g[::-1]).conj() + Dra
        gDown = np.array([np.dot(r ,a.conj())])
        return np.concatenate((gUp, gDown))

    def _Burg(self):
        """
        Uses the Standard version of Burg Algorithm to compute the power spectral
        density. The order is selected by the minimization of the chosen method

        Returns
        -------
        P: 'np.float'                  
            Variance of white noise for the associated autoregressive process
                                       
        a_k: 'np.ndarray'             
            The coefficient used to compute the power spectral density (Shape (N,)) 
            
        optimization: 'np.ndarray'    
            The values of the chosen optimisation_method at every iteration 
            (Shape (N,))   
        """
        
        #initialization of variables
        P_0 = (self.data ** 2).mean()
        P = [P_0]
        a_0 = 1
        a_k = [np.array([a_0])]
        _f = np.array(self.data)
        _b = np.array(self.data)
        optimization = []
        early_stop_step = self.mmax/50
        idx = None
        old_idx = 0
        #Burg's recursion
        for i in range(self.mmax):
            f = _f[1:]
            b = _b[:-1]
            den = np.dot(f, f) + np.dot(b, b)
            k = - 2 * np.dot(f.T, b) / den
            a_k.append(self._updatePredictionCoefficient(a_k[i], k))
            P.append(P[i] * (1 - k * k.conj()))
            _f = f + k * b
            _b = b + k * f
            #print('P: ', P, '\nak: ', a_k[-1])
            optimization.append(self._optimizer(P, a_k[-1], self.N, i + 1))
                #checking if there is a minimum (every some iterations) if early_stop option is on
            if ((i % early_stop_step == 0 and i !=0) or (i >= self.mmax-1)) and self.early_stop:
                idx = np.argmin(optimization) + 1
                if old_idx < idx:
                    old_idx = idx
                else:
                    old_idx = idx
                    break
        if not self.early_stop:
            idx = np.argmin(optimization) + 1

        return P[idx], a_k[idx], optimization
    
    def _updatePredictionCoefficient(self, x, reflectionCoefficient):
        """
        Uses the Levinson recursion to update the prediction error coefficients

        Parameters
        ----------
        x : 'np.ndarray'
            The i'th order forward prediction error coefficients (Shape (i,)).
        reflectionCoefficient : 'np.float'
            The i'th order reflection coefficients, used to update the forward
            prediction error coefficients via the solution of Levinson Recursion.

        Returns
        -------
        'np.ndarray'
            The updatet forward prediction error coefficients at order i + 1 (Shape (i + 1,)).

        """
        new_x = np.concatenate((x, np.zeros(1)))
        return new_x + reflectionCoefficient * new_x[::-1]
    
    def forecast(self, length, number_of_simulations, P = None, include_data = False): 
        """
        Forecasting on the observed process for a total number of points given 
        by length. It computes number_of_simulations realization of the forecast time series.
        This method can only be used if a_k coefficients have been computed 
        already. Use solve method before forecasting. 

        Parameters
        ----------
        length : 'np.int'
            Number of future points to be predicted 
            
        number_of_simulations : 'np.int'
            Total number of simulations of the process
            
        P : 'np.float'
            Variance of white noise for the autoregressive process. 
            Default is None and uses the estimate obtained with Burg's algorithm.
            
        include_data: `bool`
            Whether to prepend to the output the input time series

        Returns
        -------
        predictions : 'np.ndarray'
            Array containing the forecasted points for every simulation of the
            process (Shape (number_of_simulations, length))

        """
        if self.P is None or self.a_k is None:
            raise RuntimeError("PSD analysis is not performed yet: unable to forecast the data. You should call solve() before forecasting")
        if P is None: P = self.P 
        p = self.a_k.size - 1 
        predictions = np.zeros((number_of_simulations, p + length))
        predictions[:,:p] = self.data[-p:]
        coef = - self.a_k[1:][::-1]
        for i in range(length): 
            sys.stderr.write('\r {0} of {1}'.format(i + 1, length))
            predictions[:, p + i] = predictions[:, i: p + i] @ coef +\
                         np.random.normal(0, np.sqrt(P), size = number_of_simulations)
        sys.stderr.write('\n')
        if not include_data:
            return predictions[:,p:]
        return predictions

    def generate_noise(self, T, sampling_rate = 1., fmin = None, fmax = None, N_series = 1):
        """
        Generate some noise starting from a mesa object. The noise generated has the same features as the data given in input to the mesa.
            
        Parameters
        ----------
                T: `float`                  Length (in seconds) of the signal to generate
                sampling_rate: `np.float`   Sampling rate for the time series to generate
                fmin: `float`               Minimum frequency in the signal (if None, is equal to zero)
                fmax: `float`               Maximum frequency in the signal (if None, Nyquist frequency is used: f_Ny = 0.5*sampling_rate)
                N_series `float`            Number of time series to generate
            
            Returns:
            --------
                times: `np.ndarray`                  time grid at which the noise time series is evaluated at
                times_series: `np.ndarray`           time series (shape (N_series, sampling_rate*T) )
                frequencies: `np.ndarray`            frequency grid at which the noise frequency series is evaluated at
                frequency_series: `np.ndarray`       frequency series from which the times_series is computed (shape (N_series, sampling_rate*T) )
                psd: `np.ndarray`                    PSD of the model, including both positive and negative frequencies
        """
        df      = 1 / T
        N       = int(sampling_rate * T)
        times   = np.linspace(0., T , N) 
        if fmin == None: fmin = 0
        if fmax == None: fmax = (N / 2) / T
        # filter out the bad bits
        kmin = np.int(fmin/df)
        kmax = np.int(fmax/df) + 1
        
        # generate the FD noise
        frequencies = df * np.linspace(kmin, kmax, int(N / 2 + 1)) #(D,) #df * N / 2 is Ny frequency, + 1 needed because arange cuts last term
        psd = self.spectrum(1/sampling_rate, frequencies)

        sigma = np.sqrt(psd /  df * .5) #(D,)
        phi = np.random.uniform(0, 2 * np.pi, len(sigma))
        frequency_series = np.einsum('ij,j -> ij',np.random.normal(0, 1, (N_series,len(sigma))) + 1j * np.random.normal(0, 1, (N_series,len(sigma)) ), sigma) #(N_series,D)
          
        # inverse FFT to return the TD strain
        time_series = np.fft.irfft(frequency_series) * df * N ##(N_series, N )
        return times, np.squeeze(time_series), frequencies, np.squeeze(frequency_series), psd

