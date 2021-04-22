import sys
sys.path.insert(0, '..')
from memspectrum import MESA, loss_function
from GenerateTimeSeries import generate_data 
import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt

data = np.loadtxt('../../GWAnomalyDetection/maxent/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt.gz')
#data = np.loadtxt('../climate_data/data_@1h_len95689.dat')

plot_type = 'LIGO' #'climate' #'LIGO'
srate = 1. #4096.

m_FPE = MESA()
m_LL = MESA()

print(plot_type)

if False:
	P, a_k_FPE, opt_FPE = m_FPE.solve(data, optimisation_method = 'FPE', verbose = False)
	m_FPE.save('LL_validation/FPE_model_{}'.format(plot_type))
	np.savetxt('LL_validation/FPE_opt_{}'.format(plot_type), opt_FPE)

	P, a_k_LL, opt_LL = m_LL.solve(data, optimisation_method = 'LL', early_stop = False, m = 10000, verbose = True)
	m_LL.save('LL_validation/LL_model_{}'.format(plot_type))
	np.savetxt('LL_validation/LL_opt_{}'.format(plot_type), opt_LL)
else:
	m_LL.load('LL_validation/LL_model_{}'.format(plot_type))
	m_FPE.load('LL_validation/FPE_model_{}'.format(plot_type))
	a_k_FPE = m_FPE.a_k
	a_k_LL = m_LL.a_k
	opt_LL = np.loadtxt('LL_validation/LL_opt_{}'.format(plot_type))
	opt_FPE = np.loadtxt('LL_validation/FPE_opt_{}'.format(plot_type))

l = loss_function('LL')
l._set_data(data)
PSD_LL = m_LL.spectrum( 1./srate, onesided = False)[1]
PSD_FPE = m_FPE.spectrum( 1./srate, onesided = False)[1]
print("LL")
l(0,0,0,0,PSD_LL)
print("FPE")
l(0,0,0,0,PSD_FPE)

#quit()

plt.figure()
plt.loglog(*m_LL.spectrum( 1./srate, onesided = True), label = 'LL', zorder = 10)
plt.loglog(*m_FPE.spectrum( 1./srate, onesided = True), label = 'FPE',zorder = 0)
plt.legend()
plt.savefig('LL_validation/PSDs_{}.pdf'.format(plot_type))

fig, ax = plt.subplots(2,1, sharex = True)

ax[0].plot(range(len(opt_LL)), opt_LL/np.abs(opt_LL[0]), label = 'LL')
ax[0].axvline(len(a_k_LL), c = 'r')
ax[1].plot(range(len(opt_FPE)), opt_FPE/np.abs(opt_FPE[0]), label = 'FPE')
ax[1].axvline(len(a_k_FPE), c = 'r')
ax[1].set_yscale('log')
plt.legend()
plt.savefig('LL_validation/loss_values_{}.pdf'.format(plot_type))
plt.show()

