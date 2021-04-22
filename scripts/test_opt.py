import sys
sys.path.insert(0, '..')
from memspectrum import MESA 
from GenerateTimeSeries import generate_data 
import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt

#data = np.loadtxt('../../GWAnomalyDetection/maxent/data/H-H1_GWOSC_4KHZ_R1-1126259447-32.txt.gz')
data = np.loadtxt('../climate_data/data_@1h_len95689.dat')

plot_type = 'climate' #'LIGO'

m_FPE = MESA()
P, a_k_FPE, opt_FPE = m_FPE.solve(data, optimisation_method = 'FPE', verbose = False)
m_FPE.save('LL_validation/FPE_model_{}'.format(plot_type))

m_LL = MESA()
P, a_k_LL, opt_LL = m_LL.solve(data, optimisation_method = 'LL', early_stop = False, m = 60000, verbose = False)
m_LL.save('LL_validation/LL_model_{}'.format(plot_type))

plt.figure()
plt.loglog(*m_LL.spectrum( 1./4096., onesided = True), label = 'LL')
plt.loglog(*m_FPE.spectrum( 1./4096., onesided = True), label = 'FPE')
plt.legend()
plt.savefig('LL_validation/PSDs_{}.pdf'.format(plot_type))

fig, ax = plt.subplots(2,1, sharex = True)

ax[0].plot(range(len(opt_LL)), opt_LL/np.abs(min(opt_LL)), label = 'LL')
ax[0].axvline(len(a_k_LL), c = 'r')
ax[1].plot(range(len(opt_FPE)), opt_FPE/np.abs(min(opt_FPE)), label = 'FPE')
ax[1].axvline(len(a_k_FPE), c = 'r')
ax[1].set_yscale('log')
plt.legend()
plt.savefig('LL_validation/loss_values_{}.pdf'.format(plot_type))
plt.show()

