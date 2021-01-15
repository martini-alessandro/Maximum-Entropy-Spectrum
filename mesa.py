"""
mesa
====

A nice package for Maximum Entropy Spectral Analysis.
Include short description here: it will appear when typing help(mesa)

"""

import sys
import numpy as np
import warnings

#Stefano: put as everywhere as possible the dimensions of the variables we are using?
#Stefano: make the help of every function?
#Stefano: should the sampling rate be given at initialization? Or is it fine like it is today?
#Stefano: why don't we merge generate Time series in the mesa class? So we have a single file

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
        Various optimizers that can be used to infer the best recursive order
        to estimate the power spectral density

        Parameters
        ----------
        method : 'str'
            Available methods are 'FPE', 'OBD', 'CAT', 'AIC' or 'Fixed'. 

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
        Akaike method to choose the most representative recursive method. 

        Parameters
        ----------
        P : 'np.float'
            Current estimate value for the variance of the process.
        N : 'np.int'
            Length of the dataset.
        m : 'np.int'
            Current order of the recursion.

        Returns
        -------
        'Float'
            The value of Akaike's optimizer.

        """
        return P[-1] * (N + m + 1) / (N - m - 1)
    
    def _AIC(self, P, N, m):
        """
        Second Akaike method to choose the most representative recursive method. 
        It is asymptotically equivalent to "FPE"

        Parameters
        ----------
        P : 'np.float'
            Current estimate value for the variance of the process.
        N : 'np.int'
            Length of the dataset.
        m : 'np.int'
            Current order of the recursion.

        Returns
        -------
        'Float'
            The value of second Akaike's optimizer. 

        """
        return np.log(P[-1]) + (2 * m) / N
    
    def _CAT(self, P, N, m):
        """
        Parzen method to choose the most representative recursive method. 
        It is asymptotically equivalent to "FPE"

        Parameters
        ----------
        P : 'np.float'
            Current estimate value for the variance of the process.
        N : 'np.int'
            Length of the dataset.
        m : 'np.int'
            Current order of the recursion.

        Returns
        -------
        'Float'
            The value of Parzen's optimizer. 

        """
        if m == 0:
            return np.inf
        P = np.array(P[1:])
        k = np.linspace(1, m, m)
        PW_k = N / (N - k) * P
        return 1 / (N * PW_k.sum())- (1 / PW_k[-1])
    
    def _OBD(self, P, a_k, N, m):
        """
        Rao method to choose the most representative recursive method. 
        It is asymptotically equivalent to "FPE"

        Parameters
        ----------
        P : 'np.float'
            Current estimate value for the variance of the process.
        a_k: 'np.ndarray'
            Current estimate for the autoregressive coefficients. (Shape (N,))
        N : 'np.int'
            Length of the dataset.
        m : 'np.int'
            Current order of the recursion.

        Returns
        -------
        'Float'
            The value of second Akaike's optimizer. 

        """
        P_m = P[-1]
        P = np.array(P[:-1])
        return (N - m - 2)*np.log(P_m) + m*np.log(N) + np.log(P).sum() + (a_k**2).sum()
    
    def _Fixed(self, m):
        return 1./m

class MESA(object):
    """
    Class that implements Burg Method to compute the power spectral density
    
    init: data: `np.ndarray` shape (N)
    
    """
    def __init__(self, data):
        
        self.data = data
        self.N    = len(self.data)
        
    def _spectrum(self, dt, N):
        """
        Computes the Power Spectral density of the model. The PSD is evaluated on a standard grid of frequency (as given by np.fft.fftfreq), without postprocessing.
        
        Parameters:
        ----------
        dt: `np.float`      
            Sampling rate for the time series
        N: `np.float`       
            Length of the frequency grid
        
        Returns 
        ----------
        spectrum: `np.ndarray`   
            PSD of the model, including both positive and negative frequencies (shape (N,))
        freq: `np.ndarray`       
            frequencies at which spectrum is evaluated (as provided by np.fft.fftfreq) (shape (N,))
        
        """
        den = np.fft.fft(self.a_k, n=N)
        spectrum = dt * self.P / (np.abs(den) ** 2)
        return spectrum, np.fft.fftfreq(N,dt)


    def spectrum(self, dt, frequencies = None): 
        """
        Computes the power spectral density of the time series. It can be evaluated
        on standard sampling frequencies interval or a chosen interval 

        Parameters
        ----------
        dt :  'np.float'
            The sampling rate of the time series 
        frequencies : 'np.ndarray', optional
            frequencies at which power spectral density has to be computed. 
            The default is None and uses standard sampling frequencies. (shape (N,))



        Returns
        -------
        If frequency is None: 
            spectrum: 'np.ndarray'
                The Power spectral density  (Shape (N,))
            frequencies: 'np.ndarray'
                The sampling frequencies (Shape (N,))
                
        If frequency is np.ndarray
            spectrum: 
                The spectrum computed on the passed frequencies (Shape (N,))

        """

        f_ny = .5 / dt    #Nyquist frequency (minimum frequency that can be resolved with a given sampling rate 1/dt)
        
        if type(frequencies) == np.ndarray and np.max(frequencies) > f_ny + 0.1:
            warnings.warn("Some of the required frequencies are higher than the Nyquist frequency ({} Hz): a zero PSD is returned there.".format(f_ny))
        if frequencies == None: 
            N = self.N
        elif isinstance(frequencies, np.ndarray): 
            df = np.min(np.abs(np.diff(frequencies))) * 0.9 #minimum precision required by the user given grid (*0.9 to be safe)
            N = int(2. * f_ny / df)
        spec, f_spec = self._spectrum(dt, N)    
        if frequencies == None: 
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
        

        Parameters
        ----------
        m : 'np.int', optional
            The maximum number of recursions to be performed by the algorithm.
            If no argument is passed it is chosen to be m = 2N / log(2N). 
            The default is None.
            
        optimisation_method : 'str', optional
            Chose the optimizer to be minimized in order to select the best 
            recursive order. Available methods are "FPE", "CAT", "OBD", "AIC" 
            or "Fixed". If "Fixed" the order is equivalent to m. 
            The default is "FPE".
            
        method : 'str', optional
            Selects the algorithm for the computation of the power spectral density.
            Available choices are: 'Fast', 'Standard'. The default is "Fast".
            
        regularisation : 'np.float', optional
            Tikhonov regularisation, should be a number slightly larger than one
            If nothing is passed, no regularisation is implemented. The default is 1.0.
            
        early_stop : 'Boolean', optional
            When early stop, the recursion broke if no new global minimum is found
            after 200 iterations. It is reccomended with every optimiser but CAT.
            The default is True.

        Returns
        -------
        P: 'np.float'
            The variance of the white noise error of the associated AR process 
        a_k: 'np.array'
            The coefficients used to compute the power spectral density 
        Optimization: 'np.array'
            The value of the chosen optimizer for every iteration. (shape (N,))

        """
        
        self.regularisation = regularisation
        self.early_stop = early_stop
        if m == None:
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
        Implement the Fast Burg algorithm for the computation of the power spectral
        density

        Returns
        -------
        P: 'np.float'
            The variance of the white noise error of the associated AR process 
        a_k: 'np.array'
            The coefficients used to compute the power spectral density 
        Optimization: 'np.array'
            The value of the chosen optimizer for every iteration.

        """
        #FIXME: if we decide to keep the early stop option, it must be a parameter for function self.solve

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
        Implement the Fast Burg algorithm for the computation of the power spectral
        density

        Returns
        -------
        P: 'np.float'
            The variance of the white noise error of the associated AR process 
        a_k: 'np.array'
            The coefficients used to compute the power spectral density 
        Optimization: 'np.array'
            The value of the chosen optimizer for every iteration.

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
        new_x = np.concatenate((x, np.zeros(1)))
        return new_x + reflectionCoefficient * new_x[::-1]
    
    def forecast(self, length, number_of_simulations, P = None):
        #FIXME: can this be vectorized?
        if P == None: P = self.P
        p = self.a_k.size - 1 
        coef = - self.a_k[1:][::-1]
        future = [] 
        for _ in range(number_of_simulations):
            #sys.stderr.write('\r%f' %((_ + 1)/number_of_simulations))
            predictions = self.data[-p:]            
            for i in range(length):#FIXME: progress bar?
               # sys.stderr.write('\r {0} of {1}'.format(i + 1, length))
                prediction = predictions[-p:] @ coef +\
                             np.random.normal(0, np.sqrt(P))
#                while prediction < 0:
#                    prediction = predictions[-p:] @ coef +\
#                             np.random.normal(0, np.sqrt(P))
                predictions = np.append(predictions, prediction)
            future.append(predictions[p:])
        #sys.stderr.write('\n')
        return np.array(future)
    
    def forecast_vectorized(self, length, number_of_simulations, P = None, include_data = False): 
        """
        It uses the AR coefficients computed solving Burg's Algorithm to forecast
        on the observed time series. It can only be used if solve method has been used
        already. 

        Parameters
        ----------
        length : 'np.int'
            Number of points to be predicted 
            
        number_of_simulations : 'np.int'
            Total number of simulation to be performed
            
        P : 'np.float', optional
            The variance of white noise of the estiamted AR process. If None is 
            passed P is the one computed solving the method. The default is None.
            
        include_data : 'Booelan', optional
            Return the first p data used for forecasting. The default is False.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if not isinstance(P,float): P = self.P 
        p = self.a_k.size - 1 
        predictions = np.zeros((number_of_simulations, p + length))
        print(predictions.shape,p, self.P)
        predictions[:,:p] = self.data[-p:]
        coef = self.a_k[1:][::-1]
        for i in range(length): 
            sys.stderr.write('\r {0} of {1}'.format(i + 1, length))
            predictions[:, p + i] = predictions[:, i: p + i] @ coef +\
                         np.random.normal(0, np.sqrt(P), size = number_of_simulations)
        sys.stderr.write('\n')
        if include_data:
            return predictions
        else:
            return predictions[:,p:]



#FIXME: methods to save and load psd and AR coefficients in various formats, depending on the necessary application
#FIXME: save/load class as pickle
