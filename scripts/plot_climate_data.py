import numpy as np
import matplotlib.pyplot as plt
from style_sheet import init_plotting

import sys
sys.path.insert(0,'..')
import memspectrum

plot_dir = '../paper/Images/climate_plots/'

data = np.loadtxt('../climate_data/data_@1h_len271009_1989-2020.dat')

	#spectrum computation
srate = 1./3600 #in Hz 1/(1 hours)
M_CAT = memspectrum.MESA()
M_FPE = memspectrum.MESA()
M_CAT.load("../climate_data/model_CAT_long")
M_FPE.load("../climate_data/model_FPE_long")

N = len(data)
spec_FPE, f_FPE = M_FPE.spectrum(1./srate)
spec_FPE = spec_FPE[:int(N/2)]
f_FPE = f_FPE[:int(N/2)]

spec_CAT, f_CAT = M_CAT.spectrum(1./srate)
spec_CAT = spec_CAT[:int(N/2)]
f_CAT = f_CAT[:int(N/2)]

print(len(M_CAT.a_k),len(M_FPE.a_k))

	#Forecasting
N_forecast = 24*2*365

if False:
	forecast = M_CAT.forecast(data[:len(M_CAT.a_k)],N_forecast, 100, verbose = True)
	np.savetxt('forecast_climate_data.dat', forecast)
else:
	forecast = np.loadtxt('forecast_climate_data.dat')
l,m,h = np.percentile(forecast,[5,50,95],axis = 0)
true = data[len(M_CAT.a_k):len(M_CAT.a_k)+N_forecast]

	
	#plot forecasting
fig = init_plotting()
ax = fig.gca()

t_grid = np.linspace(0,N_forecast/24.,N_forecast)
ax.plot(t_grid, (true-m), c = 'b', zorder = 1)
ax.fill_between(t_grid, true-h, true-l, color = 'r',alpha = 0.3, zorder = 0)
#ax.plot(t_grid, (true-m)/(h-l), c = 'b', zorder = 1)
#ax.plot(np.linspace(0,N_forecast,N_forecast), true-h, c = 'r')
#ax.plot(np.linspace(0,N_forecast,N_forecast), true- l, c = 'r')
ax.set_xlabel("Time (days)")
ax.set_ylabel(r"$T - T_{forecast} (K)$")
plt.tight_layout()
plt.savefig(plot_dir+"forecast_accuracy.pdf")

	#plot spectrum
fig = init_plotting()
ax = fig.gca()

unit_shift = 3600*24 #frequency is made 1/day
ax.loglog(f_CAT * unit_shift , spec_CAT, label = "CAT", c = 'b', zorder = 1)
ax.loglog(f_FPE * unit_shift , spec_FPE, label = "FPE", c = 'r', zorder = 2)

ax.axvline(1./(3600.*24.*365.)*unit_shift,lw = .5, c = 'b', ls = '--', zorder = 0)
ax.axvline(1./(3600.*24.)*unit_shift,lw = .5, c = 'b', ls = '--', zorder = 0)
ax.set_xlabel("frequency (1/day)")
ax.set_ylabel("PSD (1/Hz)")

ax.legend(loc = 'upper left')


	#inset
axins =  ax.inset_axes([0.25, 0.65, 0.75, 0.35])
ids_CAT = np.where(f_CAT>.75/unit_shift)
ids_FPE = np.where(f_FPE>.75/unit_shift)
axins.loglog(f_CAT[ids_CAT] * unit_shift , spec_CAT[ids_CAT], label = "CAT", c = 'k')
axins.loglog(f_FPE[ids_FPE] * unit_shift, spec_FPE[ids_FPE], label = "FPE", c = 'r')
axins.set_yticks([])
axins.set_yticklabels([])
#axins.set_xlabel("Frequency (1/day)")
tick_list = [str(i) for i in range(1,13)]
tick_list[8] = ''
tick_list[10] = ''
axins.set_xticklabels(tick_list)
axins.set_xticks([i for i in range(1,13)])


fig.tight_layout()
plt.savefig(plot_dir+"temp_spectrum.pdf")

plt.show()

