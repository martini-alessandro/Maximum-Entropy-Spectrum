"""
All the examples in the docs. To run at once...
"""

#Intro plot - temperature timeseries
import memspectrum
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker

timeseries = np.loadtxt('data/temperature_data_@1h_len5000.dat')
dt, T = 1, len(timeseries) #Time step, Lenght of the timeseries
N_tstep = 100
f = np.arange(1/(10*24), 0.5, 1/T)
f_values = [1/(7*24), 1/24, 1/12, 1/8, 0.25, 0.5, 1]
f_names = ['1/week', '1/day', '2/day', '3/day', '1/4hours', '1/2hours','1/hour']

M = memspectrum.MESA()
M.solve(timeseries) #perform the analysis on the given time series (a real/complex np.array)
PSD = M.spectrum(dt, f) #evaluate the PSD on the given frequency grid
future_timeseries = M.forecast(timeseries[:1000], N_tstep) #forecast from the time series

plt.figure(figsize = (3.54*2, 3.54))
plt.loglog(f,PSD)
plt.xlabel('Frequency (1/hour)')
for f_val in f_values:
	plt.axvline(f_val, ls = '--', c= 'k')
ax = plt.gca()

ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator((f_values)))
ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter((f_names)))
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('../docs/img/temperature_plot_intro.png', dpi = 200)

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

## Computing the autocorrelation
import scipy.signal
R_t_empirical = scipy.signal.correlate(time_series, time_series, mode = 'same')
R_t_empirical /= np.max(R_t_empirical)
R_t_mesa = m.compute_autocorrelation(1/srate, normalize = True, scipy_convention = True)

fig = plt.figure()
ax = plt.gca()
plt.plot(np.arange(-len(R_t_empirical)//2, len(R_t_empirical)//2)*dt, R_t_empirical, label = 'empirical')
plt.plot(np.arange(-len(R_t_mesa)//2, len(R_t_mesa)//2)*dt, R_t_mesa, label = 'mesa')
plt.axvline(m.p*dt, ls = '--', c = 'k')
plt.axvline(-m.p*dt, ls = '--', c = 'k')
plt.legend()
plt.xlabel('Time lag (s)')
left, bottom, width, height = [0.18, 0.5, 0.31, 0.31]
ax_ins = fig.add_axes([left, bottom, width, height])
ax_ins.plot(range(-m.p, m.p), R_t_empirical[len(R_t_empirical)//2-m.p:len(R_t_empirical)//2+m.p])
ax_ins.plot(range(-m.p, m.p), R_t_mesa[len(R_t_mesa)//2-m.p:len(R_t_mesa)//2+m.p])
ax_ins.set_xticks([-m.p, -m.p/2, 0, m.p/2, m.p])
ax_ins.set_xticklabels(['-p', '-p/2', '0', 'p/2', 'p'])
plt.savefig('../docs/img/autocorrelation.png')
ax.set_title('Autocorrelation')

## Forecasting
id_start = 10000
id_end = id_start + 3000
length_forecasting = 1000
times_forecasting, forecast_baseline = t_grid[id_start:id_end], time_series[id_start:id_end]
forecast = m.forecast(forecast_baseline[:-length_forecasting], length = length_forecasting, number_of_simulations = 2000, include_data = False) 
median = np.median(forecast, axis = 0) #Ensemble median 
p5, p95 = np.percentile(forecast, (5, 95), axis = 0) #90% credibility boundaries

plt.figure()#figsize = [10.4, 4.8])
plt.plot(times_forecasting[:-length_forecasting], forecast_baseline[:-length_forecasting], color = 'k')
plt.fill_between(times_forecasting[-length_forecasting:], p5, p95, color = 'b', alpha = .5, label = '90% Cr.') 
plt.plot(times_forecasting[-length_forecasting:], forecast_baseline[-length_forecasting:], color = 'k', linestyle = '-.', label = 'Observed data') 
plt.plot(times_forecasting[-length_forecasting:], median, color = 'r', label = 'median estimate')
plt.xlabel('Time (s)')
plt.savefig('../docs/img/forecasting.png')
plt.title('Time series + Forecasting')


plt.show()	

