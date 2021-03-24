import numpy as np
import matplotlib.pyplot as plt
from style_sheet import init_plotting

import sys
sys.path.insert(0,'..')
import memspectrum

import h5py

plot_dir = '../paper/Images/forecast_LIGO/'

#Time line: GPS:1164603392	UTC: 2016-12-01T04:56:15

data = h5py.File('L-L1_GWOSC_O2_4KHZ_R1-1164603392-4096.hdf5', 'r')
data = np.array(data['strain']['Strain'])

	#spectrum computation
srate = 4096 #in Hz
T_train = 1000
M = memspectrum.MESA()
if False:
	M.solve(data[:int(T_train*srate)], optimisation_method = 'CAT', method = 'Fast' , verbose = True, m = 100000, early_stop = True)
	M.save('LIGO_model_CAT')
else:
	M.load('LIGO_model_CAT')

	#Forecasting
N_forecast = srate*10 #100 s of forecasting

if False:
	forecast = M.forecast(data[int(T_train*srate)-len(M.a_k):int(T_train*srate)],N_forecast, 100, verbose = True)
	np.savetxt('forecast_LIGO_data.dat', forecast)
else:
	forecast = np.loadtxt('forecast_LIGO_data.dat')
l,m,h = np.percentile(forecast,[5,50,95],axis = 0)
sigma = np.std(forecast,axis = 0)
true = data[int(T_train*srate):int(T_train*srate)+N_forecast]
#quit()

	#plot forecasting
fig = init_plotting()
plt.close()
fig, (ax,ax2) = plt.subplots(nrows=2, sharex=True)

t_grid = np.linspace(0,N_forecast/srate,N_forecast)
ax.plot(t_grid, -(true-m), c = 'b', zorder = 1)
ax.fill_between(t_grid, -true+h, -true+l, color = 'r',alpha = 0.3, zorder = 0)
#ax.plot(t_grid, (true-m)/(h-l), c = 'b', zorder = 1)
#ax.plot(np.linspace(0,N_forecast,N_forecast), true-h, c = 'r')
#ax.plot(np.linspace(0,N_forecast,N_forecast), true- l, c = 'r')
ax.set_ylabel(r"$h_{forecast} -h$")

ax2.plot(t_grid, sigma, c = 'b', zorder = 1)


ax2.set_xlabel("Time (s)")
ax2.set_ylabel(r"$\sigma$")

axins =  ax2.inset_axes([0.6, 0.15, 0.5, 0.35])
N_inset = 2500
axins.plot(t_grid[:N_inset], -(true-m)[:N_inset], c = 'b', zorder = 1)
axins.fill_between(t_grid[:N_inset], (-true+h)[:N_inset], (-true+l)[:N_inset], color = 'r',alpha = 0.3, zorder = 0)
axins.set_ylabel(r"$h_{forecast} -h$")

plt.tight_layout()
plt.savefig(plot_dir+"forecast_accuracy.pdf")

	#plot spectrum
fig = init_plotting()
ax = fig.gca()

N = M.N
spec, f = M.spectrum(1./srate)
spec = spec[:int(N/2)]
f = f[:int(N/2)]

ax.loglog(f , spec, label = "CAT", c = 'b', zorder = 1)

ax.set_xlabel("frequency (Hz)")
ax.set_ylabel("PSD (1/Hz)")


plt.show()

