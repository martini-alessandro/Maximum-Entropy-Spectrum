memspectrum
===========

**Authors** Alessandro Martini, Stefano Schmidt, Walter del Pozzo

**email** martini.alessandr@gmail.com

**Copyright** Copyright (C) 2020 Alessandro Martini

**Licence** CC BY 4.0

**Version** 1.0.0

MAXIMUM ENTROPY ESTIMATION ALGORITHM FOR FAST PSD COMPUTATION
=============================================================

``memspectrum`` is a package for the computation of power spectral densities of a time series. 
It implements a Fast verions of Burg method of Maximum Entropy Spectral Analysis.
The method is fast and reliable and provides better performance than other standard methods.
 
The computation of the power spectral density requires solving the Levinson recursion for the 
forward prediction error coefficients a_k.
The knowledge of such coefficients allows to characterize the observed process in terms of 
an autoregressive process of order p (AR(p)), being p + 1 the lenght of the a_k array. Together
with a_k coefficients a P coefficient is estimated, and is to be interpreted as the variance of 
white noise component for the process. 
The computation of these quantities allow to perform high quality forecast for the time series.
The estimate of the autoregressive order is via a second class, that implements several methods
available in literature. 

Usage of memspectrum
====================

To get the PSD computed, the following steps are required

+ Import the data
+ Call ``MESA`` class passing data as argument

::

	from memspectrum import MESA
	m = MESA(data)

+ Compute the coefficients via the ``solve()`` method: MANDATORY for further computations 

::

	m.solve()

+ At this point you can compute the spectrum and forecast

::

	m.spectrum()
	m.forecast()

Sinusoid example 
================
To compute (and plot) the spectrum of a (noisy) sinusoidal signal:
::

	from memspectrum import MESA 
	import numpy as np

Generating the data: 
::

	N, dt = 1000, .01  #Number of samples and sampling interval
	time = np.arange(0, N) * dt
	frequency = 2  
	data = np.sin(2 * np.pi * frequency * time) + np.random.normal(.4, size = 1000) 
	plt.plot(time, data, color = 'k') 
	
.. image:: https://raw.githubusercontent.com/martini-alessandro/Maximum-Entropy-Spectrum/main/memspectrum_package/ReadMeFigures/Data.jpeg
   :width: 700px
   
   
   
Solving MESA is needed to compute PSD or forecast. 
::

	M = MESA(data) 
	M.solve() 
	
The spectrum can be computed on sampling frequencies (automatically generated) or on 
some given interval 
::

	spectrum, frequencies = M.spectrum(dt)  #Computes on sampling frequencies 
	user_frequencies = np.linspace(1.5, 2.5)
	user_spectrum = M.spectrum(dt, user_frequencies) #Computes on desired window
	
Plotting the two the following is obtained: 

.. image:: https://raw.githubusercontent.com/martini-alessandro/Maximum-Entropy-Spectrum/main/memspectrum_package/ReadMeFigures/Spectrum.jpeg
   :width: 700px
   
   
   
It can also be used to perform forecasting. For example, we consider the first 900 points 
of the data and try to infer the upcoming signal. 1000 simulations of 100 points are performed.
Real observed data are compared with median estimate and 90% Credibility regions 
::

	M = MESA(data[:-100]) 
	M.solve() 
	forecast = M.forecast(length = 100, number_of_simulations = 1000, include_data = False) 
	median = np.median(forecast, axis = 0) #Ensemble median 
	p5, p95 = np.percentile(forecast, (5, 95), axis = 0) #90% credibility boundaries
	
	plt.plot(time[:-100], data[:-100], color = 'k')
	plt.fill_between(time[-100:], p5, p95, color = 'b', alpha = .5, label = '90% Cr.') 
	plt.plot(time[-100:], data[-100:], color = 'k', linestyle = '-.', label = 'Observed data') 
	plt.plot(time[-100:], median, color = 'r', label = 'median estimate') 
	 
 

The forecast result is: 

.. image:: https://raw.githubusercontent.com/martini-alessandro/Maximum-Entropy-Spectrum/main/memspectrum_package/ReadMeFigures/Forecast.jpeg
   :width: 700px


Generating data from PSD
============================
memspectrum.generateTimeSeries provides a function that construct a time-series with a user-given power 
spectral density. It can be called as 
:: 

	from memspectrum.generateTimeSerie import generate_data
	f, psd = import wanted psd and frequency array 
	time, time_series, frequency, frequency_series, psd = generate_data(f, psd, T, sampling_rate)
	
T represent the time length of the observation and sampling rate is equivalent to 1 / dt, with dt the sampling interval
 

Installation & documentation
============================
To install the package: ::

	pip install memspectrum

It requires ``numpy``.

On the GitHub repository, a number of examples are available to the interested user:

* `gwstrain.py <https://github.com/martini-alessandro/Maximum-Entropy-Spectrum/blob/main/examples/gwstrain.py>`_: computes the PSD on a piece of gravitational waves data and perform some forecasting
* `sunspots.py <https://github.com/martini-alessandro/Maximum-Entropy-Spectrum/blob/main/examples/sunspots.py>`_: using data from sunspots, it uses memspectrum to find an autoregressive process which describes them and forecast
* `sound_MESA.py <https://github.com/martini-alessandro/Maximum-Entropy-Spectrum/blob/main/examples/sound_MESA.py>`_: given an input audio (wav) file reproducing the sound of a waterfall, it computes the PSD and generate a synthetic noise, resembling the original one.

For more advanced use or for more information, please refer to the code documentation: ::

	import memspectrum
	help(memspectrum)
	help(memspectrum.<function_name>)

For full source code (and much more) see: https://github.com/martini-alessandro/Maximum-Entropy-Spectrum
