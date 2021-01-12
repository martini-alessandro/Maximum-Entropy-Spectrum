import numpy as np
from mesa.mesa import MESA
import matplotlib.pyplot as plt

if __name__ == "__main__":

    datafile = "../mesa/data/zuerich-monthly-sunspot-numbers-.tsv"
    _, data = np.genfromtxt(datafile, unpack = True)
    T = len(data)
    dt = 1.
    srate = 1./dt
    N = data.shape[0]
    f = np.fft.fftfreq(N, d=dt)
    t = np.arange(0,T,step=dt)
    M = MESA(data)
    P, ak, _ = M.solve(method = "Fast", optimisation_method = "FPE", m = int(2*N/(2*np.log(N))))
    print("p = {0}".format(len(ak)))
    PSD, _     = M.spectrum(dt,len(f))[:N//2]

    fig = plt.figure(1)
    ax  = fig.add_subplot(111)
    ax.loglog(f[:N//2], PSD[:N//2],'-k')
    ax.set_xlim(1e-5,srate/2.)
    ax.set_xlabel("frequency (1/months)")
    ax.set_ylabel("PSD (months)")
    plt.show()

