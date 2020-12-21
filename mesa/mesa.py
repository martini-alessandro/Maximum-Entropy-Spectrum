import sys
import numpy as np

class optimizer:

    def __init__(self, method):
        self.method = method

    def __call__(self, *args):
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
        return P[-1] * (N + m + 1) / (N - m - 1)
    
    def _AIC(self, P, N, m):
        return np.log(P[-1]) + (2 * m) / N
    
    def _CAT(self, P, N, m):
        if m == 0:
            return np.inf
        P = np.array(P[1:])
        k = np.linspace(1, m, m)
        PW_k = N / (N - k) * P
        return 1 / (N * PW_k.sum())- (1 / PW_k[-1])
    
    def _OBD(self, P, a_k, N, m):
        P_m = P[-1]
        P = np.array(P[:-1])
        return (N - m - 2)*np.log(P_m) + m*np.log(N) + np.log(P).sum() + (a_k**2).sum()
    
    def _Fixed(self, m):
        return m

class MESA(object):

    def __init__(self, data):
        
        self.data = data
        self.N    = len(self.data)
        
    def spectrum(self, dt, frequency):
        N = frequency.shape[0]
        den = np.fft.fft(self.a_k,n=N)
        spec = dt * self.P / (np.abs(den) ** 2)
        return spec

    def spectrum_bis(self, frequencies, dt):
        "Alternative to spectrum: takes in input a positive frequency grid and a dt and it evaluates the PSD on that grid."
            # df = 1/(dt*N) --> N = 2 f_ny/df
            #see https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html#numpy.fft.fftfreq
        df = np.min(np.abs(np.diff(frequencies))) *0.9 #minimum precision required by the user given grid (*0.9 to be safe)
        f_ny = 0.5/dt #Nyquist frequency (minimum frequency that can be resolved with a given sampling rate 1/dt)
        if np.max(frequencies) > f_ny:
            #here we could also raise a warning and set a 0 PSD
            raise ValueError("Some of the required frequencies are higher than the Nyquist frequency: unable to continue")
        N = int(2.*f_ny/df)

            #doing the actual computation
        den = np.fft.fft(self.a_k,n=N)
        spec = dt * self.P / np.multiply(den, den.conjugate())
        f_spec = np.fft.fftfreq(N,dt) #grid at which spec is evaluated
            #only real part of spec is used (PSD is a real function!)
            #only positive frequencies of the spectrum are computed
        f_interp = np.interp(frequencies, f_spec[:int(N/2+0.5)], spec.real[:int(N/2+0.5)], left = 0., right = 0.) 

        return f_interp


    def solve(self,
              m = None,
              optimisation_method = "FPE",
              method              = "Fast",
              regularisation      = 1.0):
        
        self.regularisation = regularisation
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
    
    def _FastBurg(self):
        #Define autocorrelation
        c = np.zeros(self.mmax + 2)
        for j in range(self.mmax + 1):
            c[j] = self.data[: self.N - j] @ self.data[j : ]
        c[0] *= self.regularisation
        #Initialize variables
        a = [np.array([1])]
        P = [c[0] / self.N]
        r = 2 * c[1]
        g = np.array([2 * c[0] - self.data[0] * self.data[0].conj() - self.data[-1] * self.data[-1].conj(), r])
        #Initialize lists to save arrays
        optimization = np.zeros(self.mmax)
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
            optimization[i] = self._optimizer(P, a[-1], self.N, i + 1)
        if self._optimizer.method == "Fixed":
            idx = self.mmax
        else:
            idx = optimization.argmin() + 1
        return P[idx], a[idx], optimization
        
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
        #initialization of variables
        P_0 = (self.data ** 2).mean()
        P = [P_0]
        a_0 = 1
        a_k = [np.array([a_0])]
        _f = np.array(self.data)
        _b = np.array(self.data)
        optimization = np.zeros(self.mmax)
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
            optimization[i] = self._optimizer(P, a_k[-1], self.N, i + 1)
        #selecting the minimum for the optimizer and recording its position
        if self._optimizer.method == "Fixed":
            idx = self.mmax
        else:
            idx = optimization.argmin()+1
        return P[idx], a_k[idx], optimization
    
    def _updatePredictionCoefficient(self, x, reflectionCoefficient):
        new_x = np.concatenate((x, np.zeros(1)))
        return new_x + reflectionCoefficient * new_x[::-1]
    
    def forecast(self, length, number_of_simulations, P = None): 
        if P == None: P = self.P
        p = self.a_k.size - 1 
        coef = - self.a_k[1:][::-1]
        future = [] 
        for _ in range(number_of_simulations):
            sys.stderr.write('\r%f' %((_ + 1)/number_of_simulations))
            predictions = self.data[-p:]            
            for i in range(length):
                sys.stderr.write('\r {0} of {1}'.format(i + 1, length))
                prediction = predictions[-p:] @ coef +\
                             np.random.normal(0, np.sqrt(P))
#                while prediction < 0:
#                    prediction = predictions[-p:] @ coef +\
#                             np.random.normal(0, np.sqrt(P))
                predictions = np.append(predictions, prediction)
            future.append(predictions[p:])
        sys.stderr.write('\n')
        return np.array(future)

def autocorrelation(x, norm = 'N'):
    N = len(x)
    X=np.fft.fft(x-x.mean())
    # We take the real part just to convert the complex output of fft to a real numpy float. The imaginary part if already 0 when coming out of the fft.
    R = np.real(np.fft.ifft(X*X.conj()))
    # Divide by an additional factor of 1/N since we are taking two fft and one ifft without unitary normalization, see: https://docs.scipy.org/doc/numpy/reference/routines.fft.html#module-numpy.fft
    if norm == 'N':
        return R/N
    elif norm == None: 
        return R
    else:
        raise ValueError('this normalization is not available')
