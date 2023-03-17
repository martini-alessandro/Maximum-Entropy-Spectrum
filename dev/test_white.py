import numpy as np
from memspectrum import MESA
import matplotlib.pyplot as plt

#######################

data = np.loadtxt('climate_data/data_@1h_len95689_2009-2020.dat')
m = MESA()
m.solve(data, method = 'standard')

white_data = m.whiten(data)

m_white = MESA()
m_white.solve(white_data, method = 'standard')

m_white_trim = MESA()
m_white_trim.solve(white_data[m.get_p():-m.get_p()], method = 'standard')

plt.loglog(*m_white.spectrum(), label = 'no trim')
plt.loglog(*m_white_trim.spectrum(), label = 'trim')
plt.legend()
plt.savefig('/home/stefano/Pictures/psd_comparison.png')
plt.show()

