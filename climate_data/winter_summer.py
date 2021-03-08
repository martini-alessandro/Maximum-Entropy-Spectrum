import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'..')
import memspectrum
import scipy.signal

sys.path.insert(0,'../scripts')
import welch


data_summer = np.loadtxt('data_@1h_len95689.dat')[(195-30)*24:(195+30)*24] #summer
data_winter = np.loadtxt('data_@1h_len95689.dat')[(2567-30)*24:(2567+30)*24] #winter

srate = 1./3600 # 1/(1 hours)
M_summer = memspectrum.MESA()
M_winter = memspectrum.MESA()
M_summer.solve(data_summer, method ='fast', optimisation_method = 'FPE', early_stop = False)
M_winter.solve(data_winter, method ='fast', optimisation_method = 'FPE', early_stop = False)

spec_summer, f = M_summer.spectrum(1./srate)
spec_winter, f = M_winter.spectrum(1./srate)
N = len(f)
spec_summer = spec_summer[:int(N/2)]
spec_winter = spec_winter[:int(N/2)]
f = f[:int(N/2)]

unit_shift = 3600*24 #frequency is made 1/day
plt.loglog(f*unit_shift,spec_summer, label = 'summer') 
plt.loglog(f*unit_shift,spec_winter, label = 'winter')
plt.legend()
plt.show()
