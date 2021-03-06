memspectrum
===========

**Authors** Alessandro Martini, Stefano Schmidt, Walter del Pozzo

**emails** martini.alessandr@gmail.com, stefanoschmidt1995@gmail.com, walter.delpozzo@ligo.org

**Copyright** Copyright (C) 2020 Alessandro Martini

**Licence** CC BY 4.0

**Version** 1.1.0.post

MAXIMUM ENTROPY ESTIMATION ALGORITHM FOR ACCURATE PSD COMPUTATION
=================================================================

``memspectrum`` is a package for the computation of power spectral densitiy (PSD) of time series. 
It implements a fast numpy verion of the Burg method for Maximum Entropy Spectral Analysis.
The method is fast and reliable and shows better performance than other standard methods.

The method is based on the maximum entropy principle, and it allows to make minimal
 assumptions on unavailable information. Furthermore, it provides a beatiful link between spectral 
 analysis and the theory of autoregressive processes

The PSD is expressed in terms of a set of coefficients a_k plus an overall scale factor P.
The a_ks are obtained recursively through the Levinson recursion.
The knowledge of such coefficients allows to characterize the observed time series in terms of 
an autoregressive process of order p (AR(p)), being p + 1 the lenght of the a_k array.
The a_k coefficients are the autoregressive coefficients, while the P scale factor can be interpreted 
as the variance of white noise component for the process. 
Once the link with an AR(p) process is established, high quality forecast for the time series is straightforward.

Usage of memspectrum
====================

To get the PSD computed, the following steps are required

+ Import the data
+ Import `memspectrum` and create an instance of ``MESA`` class:

::

	from memspectrum import MESA
	m = MESA()

+ Compute the autoregressive coefficients via the ``solve()`` method (*required* for further computations)

::

	m.solve(data)

+ At this point you can compute the spectrum and forecast N future observations

::

	spec, frequencies = m.spectrum(dt)
	predicted_data = m.forecast(data, N)

Sinusoid example 
================
To compute (and plot) the spectrum of a (noisy) sinusoidal signal:
::

	from memspectrum import MESA 
	import numpy as np
	import matplotlib.pyplot as plt

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

	M = MESA() 
	M.solve(data) 
	
The spectrum can be computed on sampling frequencies (automatically generated) or on 
some given interval 
::

	spectrum, frequencies = M.spectrum(dt)  #Computes on sampling frequencies 
	user_frequencies = np.linspace(1.5, 2.5)
	user_spectrum = M.spectrum(dt, user_frequencies) #Computes on desired frequency grid
	
The two spectra look like

.. image:: https://raw.githubusercontent.com/martini-alessandro/Maximum-Entropy-Spectrum/main/memspectrum_package/ReadMeFigures/Spectrum.jpeg
   :width: 700px
   
   
It can also be used to perform forecasting. For example, we consider the first 900 points 
of the data and try to infer the upcoming signal. 1000 simulations of 100 points are performed.
Real observed data are compared with median estimate and 90% Credibility regions 
::

	M = MESA() 
	M.solve(data[:-100]) 
	forecast = M.forecast(data[:-100], length = 100, number_of_simulations = 1000, include_data = False) 
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
Module ``memspectrum.GenerateTimeSeries`` provides a function that construct a time-series with a user-given power spectral density. It can be called as 
:: 

	from memspectrum.GenerateTimeSeries import generate_data
	f, psd = (whathever psd and frequency array you like)
	time, time_series, frequency, frequency_series, psd = generate_data(f, psd, T, sampling_rate)
	
where T represents the time length of the observation and the sampling rate is equivalent to the inverse of the sampling interval
 

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

If you feel that you need to know more about the code, or you just want to say hi, feel free to contact one of the authors.
