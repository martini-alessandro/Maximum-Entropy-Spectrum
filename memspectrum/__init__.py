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

.. code-block:: python

	import memspectrum

	M = memspectrum.MESA()
	M.solve(time_series) #perform the analysis on the given time series (a real/complex np.array)
	M.spectrum(dt,f) #evaluate the PSD on the given frequency grid
	M.forecast(data, N_tstep) #forecast from the time series

With the proper input, this will produce this nice plot:

.. image:: img/temperature_plot_intro.png

If you're curious, this is the Power Spectral Density of the historical series of temperature measured around the city of Milan (Italy) with 1 hour rate. Data are taken from `Copernicus <https://cds.climate.copernicus.eu/#!/home>`_. You can clearly see that the temperatures are very nicely correlated with themselfs on a day timescale. This amounts to a large peak at a frequency of 1 per hour and its multiples.

You can find this and many other examples around this documentation. For some example code, you can check the `example <https://github.com/martini-alessandro/Maximum-Entropy-Spectrum/tree/main/examples>`_ folder of this repo.


"""

from .memspectrum import MESA
