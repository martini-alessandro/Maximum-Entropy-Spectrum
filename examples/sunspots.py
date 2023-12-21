"""
Script to compute the power spectral density of the measured timeseries of the number of sunspots measured.
"""

from memspectrum import MESA

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    datafile = "data/zuerich-monthly-sunspot-numbers-.tsv"
    _, data = np.genfromtxt(datafile, unpack = True)
    T = len(data)
    #print([d for d in data])
    dt = 1.
    srate = 1./dt
    N = data.shape[0]
    f = np.fft.fftfreq(N, d=dt)
    t = np.arange(0,T,step=dt)
    M = MESA()
    P, ak, _ = M.solve(data, method = "Fast", optimisation_method = "FPE", m = int(2*N/(2*np.log(N))))
    print("Length of the selected filter: p = {0}".format(len(ak)))
    _, PSD     = M.spectrum(dt)[:N//2]

    fig = plt.figure(1)
    ax  = fig.add_subplot(111)
    ax.loglog(f[:N//2], PSD[:N//2],'-k')
    ax.set_xlim(1e-5,srate/2.)
    ax.set_xlabel("frequency (1/months)")
    ax.set_ylabel("PSD (months)")
    plt.show()

