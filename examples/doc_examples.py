"""
All the examples in the docs...
"""
#Getting the data
from memspectrum.GenerateTimeSeries import generate_data
import numpy as np

T, srate = 10., 4096.	
f_grid_analytic = np.linspace(.1, 4096., 1000)
psd_analytic = 1e-10*np.exp((np.log(f_grid_analytic)- np.log(20))**2)
	
t_grid, time_series, f_grid, frequency_series, interpolated_psd = generate_data(f_grid_analytic, psd_analytic, T, srate, seed = 0)

## Solving the spectrum
import memspectrum
import matplotlib.pyplot as plt

m = memspectrum.MESA()
m.solve(time_series, method = 'Standard', optimisation_method = 'VM')
f_mesa, psd_mesa = m.spectrum(dt = 1./srate, onesided = True)
	
plt.figure()
plt.plot(t_grid, time_series)
plt.xlabel('time (s)')
plt.savefig('../docs/img/timeseries_analytical.png')
plt.title('Time series')
	
plt.figure()
plt.loglog(f_grid_analytic, psd_analytic, ls = '--', label = 'true PSD')
plt.loglog(f_mesa, psd_mesa, label = 'estimated PSD')
plt.xlabel('f (Hz)')
plt.ylabel('PSD (1/Hz)')
plt.legend()
plt.savefig('../docs/img/psd_true_vs_analytical.png')
plt.title('Estimated PSD (vs true PSD)')


## Whitening

white_time_series = m.whiten(time_series, trim = 0)

plt.figure()
plt.plot(t_grid, white_time_series)
plt.xlabel('time (s)')
plt.savefig('../docs/img/white_timeseries_notrim.png')
plt.title('Whitened timeseries (no trim)')

white_time_series = m.whiten(time_series)
plt.figure()
plt.plot(t_grid[m.p:-m.p], white_time_series)
plt.xlabel('time (s)')
plt.savefig('../docs/img/white_timeseries.png')
plt.title('Whitened timeseries')

m_white = memspectrum.MESA()
m_white.solve(white_time_series)
f, PSD = m_white.spectrum(dt = 1./srate, onesided = True)
print("AR order for the whitened data: ", m_white.p)

## Generating noise
_, fake_time_series, _, _, _ = m.generate_noise(T, sampling_rate = srate)

m = memspectrum.MESA()
m.solve(fake_time_series, method = 'Standard', optimisation_method = 'VM')
f_mesa_generated, psd_mesa_generated = m.spectrum(dt = 1./srate, onesided = True)

plt.figure()
plt.loglog(f_grid_analytic, psd_analytic, ls = '--', label = 'true PSD')
plt.loglog(f_mesa, psd_mesa, label = 'data PSD')
plt.loglog(f_mesa_generated, psd_mesa_generated, label = 'fake data PSD')
plt.xlabel('f (Hz)')
plt.ylabel('PSD (1/Hz)')
plt.legend()
plt.savefig('../docs/img/psd_comparison_generate_noise.png')
plt.title('PSD estimated from syntetic noise')



## Forecasting
id_start = 10000
id_end = id_start + 3000
length_forecasting = 1000
times_forecasting, forecast_baseline = t_grid[id_start:id_end], time_series[id_start:id_end]
forecast = m.forecast(forecast_baseline[:-length_forecasting], length = length_forecasting, number_of_simulations = 2000, include_data = False) 
median = np.median(forecast, axis = 0) #Ensemble median 
p5, p95 = np.percentile(forecast, (5, 95), axis = 0) #90% credibility boundaries

plt.figure(figsize = [10.4, 4.8])
plt.plot(times_forecasting[:-length_forecasting], forecast_baseline[:-length_forecasting], color = 'k')
plt.fill_between(times_forecasting[-length_forecasting:], p5, p95, color = 'b', alpha = .5, label = '90% Cr.') 
plt.plot(times_forecasting[-length_forecasting:], forecast_baseline[-length_forecasting:], color = 'k', linestyle = '-.', label = 'Observed data') 
plt.plot(times_forecasting[-length_forecasting:], median, color = 'r', label = 'median estimate') 
plt.savefig('../docs/img/forecasting.png')
plt.title('Time series + Forecasting')


plt.show()	

