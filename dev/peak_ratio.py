import numpy as np
from memspectrum import MESA
import matplotlib.pyplot as plt
import scipy.signal

n = 9996 # number of data points
T = 1      # total duration of the data in second
fs = n/T  # sample rate
# Defining the sin wave
f1 = 1000
f2 = 300
a1 = 100
a2 = 10

# Generating some signal.
t = np.linspace(0,1,n)
s1 = a1 * np.sin(2 * np.pi * f1 * t) + a2 * np.sin(2 * np.pi * f2 * t)  + np.random.normal(0, scale=1, size=n)

m = MESA()
m.solve(s1, method = 'Standard', optimisation_method = 'VM')

psd_welch = scipy.signal.welch(s1, fs, window='tukey')


# Plotting

plt.figure()
plt.loglog(*m.spectrum(1/fs, onesided = True), label = 'mesa')
plt.loglog(*psd_welch, label = 'welch')
plt.xlabel(r"$f \;\;(Hz)$")
plt.ylabel(r"$PSD \;\;(1/\sqrt{Hz})$")
plt.legend()
plt.show()
