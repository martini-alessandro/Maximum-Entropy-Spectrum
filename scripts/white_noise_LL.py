import sys
sys.path.insert(0, '..')
from memspectrum import MESA, loss_function
from GenerateTimeSeries import generate_data 
import numpy as np
import scipy.stats 
import matplotlib.pyplot as plt

data = np.random.normal(0,1,10000)
srate = 1.

m_FPE = MESA()
m_LL = MESA()

P, a_k_FPE, opt_FPE = m_FPE.solve(data, optimisation_method = 'FPE', early_stop = False, verbose = True)

P, a_k_LL, opt_LL = m_LL.solve(data, optimisation_method = 'LL', early_stop = False, m = None, verbose = True)


plt.figure()
plt.loglog(*m_LL.spectrum( 1./srate, onesided = True), label = 'LL', zorder = 10)
plt.loglog(*m_FPE.spectrum( 1./srate, onesided = True), label = 'FPE',zorder = 0)
plt.legend()

fig, ax = plt.subplots(2,1, sharex = True)

ax[0].plot(range(len(opt_LL)), opt_LL/np.abs(opt_LL[0]), label = 'LL')
ax[0].axvline(len(a_k_LL), c = 'r')
ax[1].plot(range(len(opt_FPE)), opt_FPE/np.abs(opt_FPE[0]), label = 'FPE')
ax[1].axvline(len(a_k_FPE), c = 'r')
ax[1].set_yscale('log')
plt.legend()
plt.show()
