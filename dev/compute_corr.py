import memspectrum
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import scipy.signal

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

## Computing the autocorrelation
R_t_empirical = scipy.signal.correlate(time_series, time_series, mode = 'same')
R_t_mesa = m.compute_autocorrelation(1/srate, normalize = True)
N = len(R_t_mesa)
#R_t_mesa = np.concatenate([R_t_mesa[N//2+1:], R_t_mesa[:N//2]])


R_t_empirical /= np.max(R_t_empirical)
R_t_mesa /= np.max(R_t_mesa)

plt.figure()
plt.plot(R_t_empirical, label = 'empirical')
plt.plot(R_t_mesa, label = 'mesa')
plt.axvline(N//2+m.p, ls = '--', c = 'k')
plt.axvline(N//2-m.p, ls = '--', c = 'k')
plt.legend()
plt.show()
quit()
