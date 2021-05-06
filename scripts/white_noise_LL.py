import sys
sys.path.insert(0, '..')
from memspectrum import MESA, loss_function
from GenerateTimeSeries import generate_data 
import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt
import welch
from matplotlib.lines import Line2D


data = np.random.normal(0,1,20000)
srate = 1.

m_FPE = MESA()
m_LL = MESA()

P, a_k_FPE, opt_FPE = m_FPE.solve(data, optimisation_method = 'FPE', early_stop = False, verbose = True)

P, a_k_LL, opt_LL = m_LL.solve(data[:100], optimisation_method = 'LL', early_stop = False, m = None, verbose = True)

N_windows_list = [2, 5, 10, 100, 1000, 2000, None]
factor = 10**len(N_windows_list)

plt.figure(figsize = (6,4))
ax = plt.axes(frameon = False)
for i, N in enumerate(N_windows_list):
	if N is None:
		f_welch, PSD_welch = m_FPE.spectrum( 1./srate, onesided = True)
		ax.loglog(f_welch, PSD_welch*factor, '--',c = 'k', lw = 2)
	else:
		f_welch, PSD_welch = welch.psd(data,
		    srate,
		    srate*len(data)/N,
		    window_function  = None,
		    overlap_fraction = 0.8,
		    nfft = None,
		    return_onesided = True)
		ax.loglog(f_welch, PSD_welch/np.sqrt(2)*factor, '-',c = 'k', lw = 1)
	factor /= 10
	if i< 10:
		factor /= 10
ax.axes.get_yaxis().set_visible(False)
ax.axes.get_xaxis().set_visible(True)
xmin, xmax = ax.get_xaxis().get_view_interval()
ymin, ymax = ax.get_yaxis().get_view_interval()
ax.set_xlabel(r'$f(a.u.)$')
ax.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
plt.tight_layout()
plt.savefig('/home/stefano/Dropbox/Stefano/PhD/presentations/memspectrum_RD_call/img/welch_tuning.pdf', transparent = True)
plt.show()

plt.figure()
#plt.loglog(*m_LL.spectrum( 1./srate, onesided = True), label = 'LL', zorder = 10)
plt.loglog(*m_FPE.spectrum( 1./srate, onesided = True), label = 'FPE',zorder = 0)
plt.loglog(f_welch, PSD_welch/np.sqrt(2), '-',c = 'k')

fig, ax = plt.subplots(2,1, sharex = True)

ax[0].plot(range(len(opt_LL)), opt_LL/np.abs(opt_LL[0]), label = 'LL')
ax[0].axvline(len(a_k_LL), c = 'r')
ax[1].plot(range(len(opt_FPE)), opt_FPE/np.abs(opt_FPE[0]), label = 'FPE')
ax[1].axvline(len(a_k_FPE), c = 'r')
ax[1].set_yscale('log')
plt.legend()
plt.show()
