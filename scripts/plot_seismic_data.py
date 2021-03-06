import numpy as np
import matplotlib.pyplot as plt
from style_sheet import init_plotting

import sys
sys.path.insert(0,'..')
import memspectrum

from obspy import read

#http://rdsa.knmi.nl/fdsnws/dataselect/1/builder

#[y] = v

plot_dir = '../paper/Images/seismic_plots/'

st = read('knmi-fdsnws.mseed')

print(st)
print(st[1].stats)
print(st[1].data)

data = np.array(st[1].data)[:1000000]
data = data/1e4 #pay attention to this scaling
print("Len data (days): ", len(data)/40./(3600.*24))

del st

###################

	#spectrum computation
srate = 40. #in Hz 1/(1 hours)
M_CAT = memspectrum.MESA()
M_FPE = memspectrum.MESA()

T_train = 10000

if False:
	M_CAT.solve(data[-int(T_train*srate):], optimisation_method = 'CAT', method = 'Standard' , verbose = True, m = 10000, early_stop = True)	
	M_CAT.save('seismic_model_CAT')
	M_FPE.solve(data[-int(T_train*srate):], optimisation_method = 'FPE', method = 'Standard' , verbose = True, m = 10000, early_stop = True)	
	M_FPE.save('seismic_model_FPE')
else:
	M_CAT.load("seismic_model_CAT")
	M_FPE.load("seismic_model_FPE")

N = int(T_train*srate)
f_FPE, spec_FPE = M_FPE.spectrum(1./srate)
spec_FPE = spec_FPE[:int(N/2)]*1e8
f_FPE = f_FPE[:int(N/2)]

f_CAT, spec_CAT = M_CAT.spectrum(1./srate)
spec_CAT = spec_CAT[:int(N/2)]*1e8
f_CAT = f_CAT[:int(N/2)]

print(len(M_CAT.a_k),len(M_FPE.a_k))

	#Forecasting
N_forecast = int(40.*60*10)

if False:
	forecast = M_CAT.forecast(data[int(T_train*srate) - len(M_CAT.a_k):int(T_train*srate)], N_forecast, 100, verbose = True)
	np.savetxt('forecast_seismic_data.dat', forecast)
else:
	forecast = np.loadtxt('forecast_seismic_data.dat')
l,m,h = np.percentile(forecast[:,:N_forecast],[5,50,95],axis = 0)
sigma = np.std(forecast[:,:N_forecast], axis = 0)
true = data[int(T_train*srate):int(T_train*srate)+N_forecast]

	#plot data
fig = init_plotting()
plt.plot(np.linspace(0, len(data)/srate/60, len(data)),data*1e4)

	
	#plot forecasting
fig = init_plotting()
ax = fig.gca()

t_grid = np.linspace(0,N_forecast/srate,N_forecast)
ax.plot(t_grid, -(true-m), c = 'b', zorder = 1)
ax.fill_between(t_grid, -true+h, -true+l, color = 'r',alpha = 0.3, zorder = 0)

#ax.plot(t_grid, sigma)

ax.set_xlabel("Time (s)")
ax.set_ylabel(r"$?_{forecast} -? (??)$")
plt.tight_layout()
plt.savefig(plot_dir+"forecast_accuracy.pdf")

	#plot spectrum
fig = init_plotting()
ax = fig.gca()

ax.loglog(f_CAT  , spec_CAT, label = "CAT", c = 'b', zorder = 1)
ax.loglog(f_FPE , spec_FPE, label = "FPE", c = 'r', zorder = 2)

ax.set_xlabel("frequency (Hz)")
ax.set_ylabel("PSD (1/Hz)")

ax.legend(loc = 'upper left')

fig.tight_layout()
plt.savefig(plot_dir+"seismic_spectrum.pdf")

plt.show()

