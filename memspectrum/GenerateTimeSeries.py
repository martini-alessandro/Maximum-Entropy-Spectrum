# -*- coding: utf-8 -*-
#"""
#Module that generates a random time series with a given power spectral density
#"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def generate_data(f,
				  psd,
				  T, 
				  sampling_rate = 1.,
				  fmin = None,
				  fmax = None,
				  asd = False,
				  seed = None):
	"""
	Generate a time series with a given power spectral density.
	
	.. code::

		from memspectrum.GenerateTimeSeries import generate_data

	Parameters
	----------
	f : :class:`~numpy:numpy.ndarray`
		The frequencies over which the power spectral density is evaluated. (Shape (N,))
	psd : :class:`~numpy:numpy.ndarray`
		The power spectral density (Shape (N,))
	T : float
		The total time of the observation 
	sampling_rate : float
		The sampling rate of the output series. The default is 1..
	fmin : float
		The minimum frequency available. The default is None.
	fmax : float
		Tha maximum frequency available. The default is None.
	asd : bool
		If True, takes the square of the input power spectral density. The default is False
	seed: int
		If given, it sets a seed for the ranodom noise generation for reproducibility.

	Returns
	-------
	times: :class:`~numpy:numpy.ndarray` 
		The sampling time vector. (Shape (N,))
	time_series: :class:`~numpy:numpy.ndarray`
		The output time series  (Shape (N,))
	frequencies: :class:`~numpy:numpy.ndarray`
		The sampling frequencies (Shape (N,))
	frequency_series: :class:`~numpy:numpy.ndarray`
		The output series in frequency domain (Shape (N,))
	psd: :class:`~numpy:numpy.ndarray`
		The frequencies interpolated power spectral density (Shape (N,))

	"""
	if isinstance(seed, int): np.random.seed(seed)
	# f, psd = np.loadtxt(psd_file, unpack=True)
	if asd is True : psd = np.square(psd)
	# generate an interpolant for the PSD
	psd_int = interp1d(f, psd, bounds_error=False, fill_value='extrapolate')
	df	  = 1 / T
	N	   = int(sampling_rate * T)
	times   = np.linspace(0, T, N) 
	if fmin == None: fmin = 0
	if fmax == None: fmax = (N / 2) / T
	# filter out the bad bits
	kmin = int(fmin/df)
	kmax = int(fmax/df) + 1
	

	# generate the FD noise
	frequencies = df * np.arange(kmin, kmax) #df * N / 2 is Ny frequency, + 1 needed because arange cuts last term
	frequency_series = np.zeros(len(frequencies), dtype = np.complex128)


	sigma = np.sqrt(psd_int(frequencies) /  df * .5) 
	frequency_series = sigma * (np.random.normal(0, 1, len(sigma)) + 1j * np.random.normal(0, 1, len(sigma)))
	  
	# inverse FFT to return the TD strain
	time_series = np.fft.irfft(frequency_series, n=N) * df * N
	return times, time_series, frequencies, frequency_series, psd_int(frequencies)


