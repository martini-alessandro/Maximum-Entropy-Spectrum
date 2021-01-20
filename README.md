# Maximum Entropy Spectrum 

This repository contains an implementation of Burg's Algorithm for the computation of the power spectral density of time-series via Maximum Entropy Principle, known as Maximum Entropy Spectral Analysis (MESA). Two different methods implement the standard Algorithm and a Faster version, called FastBurg. 
The problem admits an analytical solution that can be expressed in closed form as a Discrete 
Fourier Transform of some coefficients a_k, known as forward prediction error coefficients.
The a_k coefficients are computed recursively by the use of Levinson Recursion. The role of 
MESA class is to solve the recursion. 
Various method exist to estimate the recursive order that better approximate the power spectral 
density.  A second class is inserted to implement different optimzers to choose the recursive order. 

The a_k coefficients are found to be the "best linear predictor" for the time series under study,
their computation via the former method is equivalent to a least square fitting with an autoregressive
process of order p (AR(p)). Levinson recursion also return a "P" coefficient that is equivalent to the 
variance of the white noise component for the process. This description is stationary by construction 
of Burg's Method. 

Given the a_k coefficient, they can be used to perform high quality forecasting for the future
values of the time series. A method to perform forecast is implemented in mesa class.

The code is also released as a public package `mespectrum`, available on the PyPI repository.

_________
# Code overview
## class: MESA 

	Implements the mesa algorithm, input: The time-series under study 

Call it as 
```Python
from mespectrum import MESA
M = MESA(data) 
```
#### methods: 

	`mesa.solve()`
	Solve the Levinson recursion for the computation of the a_k and P
	coefficients. 

	`mesa.spectrum()` 
	Computes the value of the power spectral density. It can be returned as both a function of the sampling frequencies or on a user-given 
	frequency grid. 
	Since P and a_k are needed, it is necessary to call .solve() method first. 
	
	`mesa.forecast()`
	It uses the forward prediction error coefficients to predict the future values of 	the time series. User can choose both the 'time len-gth' for the forecasting and the total number of simulations. 
	Since P and a_k are needed, it is necessary to call .solve() method first. 


_________

## class: optimizer
	
	Implements one of the various method to estimate the best recursive order. Once an optimizer is selected, the order is estimated as the order that minimizes the given functional relation. 
	This class is not to be called and cooperates with mesa. 
	
	

_________

## mesa.generateTimeSeries:

    Function that generate a noise-series with given power spectral 
    density

____

# A simple

A example is a sinusoidal signal with known frequency with white noise superimposition. 


#### Generating array of data:  
```Python
N, dt = 1000, .01 		#the number of samples and sampling rate 
f = 2 	 	   	    	#The frequency of the sinusoidal signal 
time = np.arange(0, N) * dt 	#The time array 
data = np.sin(2 * np.pi * f * time) + np.random.normal(scale = 0.4,
                                                       size = 1000)
```

#### Solving the recursion: 
```Python
M = mesa.mesa(data) 		#Initialize MESA class
P, a_k, opt = M.solve() 	#Solve method and returns the 
                            #coefficients and values of the opimizer
```
	
Defaul argument for solve uses Akaike's Final Prediction Error (FPE) optimizer and the fast algorithm 
	

#### Computing spectrum - two different ways 
```Python
#Computation of the spectrum on sampling frequencies 
spectrum, frequencies = M.spectrum(dt) 
	
#Computation of the spectrum on a user-given frequency grid: 
user_frequencies = np.linspace(1.5, 2.5, 1000) 
user_spectrum = M.spectrum(dt, user_frequencies) 
```


#### Forecasting and comparing result with real data 
```Python
M = MESA(data[:-100]) 
M.solve() 
forecast = M.forecast(length = 100, number_of_simulations = 1000)
```

# Installation and documentation

To install the code you might want to clone the repository

```
git clone https://github.com/martini-alessandro/Maximum-Entropy-Spectrum.git
```

Once you have it on local, you can import module mespectrum.py from your code and you're done.

You can also run some [examples](https://github.com/martini-alessandro/Maximum-Entropy-Spectrum/tree/main/examples) which provide some interesting applications of the algorithm.

The code is also released as a PyPI package:

```
pip install mespectrum
```

We did our best to document every function. If you need it, you can use the python help utility:

```
help(mespectrum.MESA.function_name)
```

For further information, feel free to email: ...
	

## References 
[J.P. Burg - Maximum Entropy Spectral Analysis](http://sepwww.stanford.edu/data/media/public/oldreports/sep06/)

[V. Fastubrg - A Fast Implementation of Burg Method](
https://svn.xiph.org/websites/opus-codec.org/docs/vos_fastburg.pdf)

## License 
[GNU General Public License v3.0](https://github.com/martini-alessandro/Maximum-Entropy-Spectrum/blob/main/LICENSE)
