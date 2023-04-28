"""
memspectrum
===========

Package that uses maximum entropy spectral Analysis to compute the spectrum 
of a given time-series. The main object is MESA, that is meant to implement Burg
Algorithm for the computation of the power spectral density.

For more information, take a look at:

	- `repository <https://github.com/martini-alessandro/Maximum-Entropy-Spectrum>`_
	- `pypy distribution <https://pypi.org/project/memspectrum/>`_
	- `documentation <https://maximum-entropy-spectrum.readthedocs.io/en/latest/>`_
	- `paper <https://arxiv.org/abs/2106.09499>`_

Basic usage (to compute the spectrum of a given time series):

::

	import memspectrum

	M = memspectrum.MESA()
	M.solve(time_series) #perform the analysis on the given time series (a real/complex np.array)
	M.spectrum(dt,f) #evaluate the PSD on the given frequency grid
	M.forecast(data, N_tstep) #forecast from the time series

"""

from .memspectrum import MESA
