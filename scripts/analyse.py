try:
    import sys
    sys.path.insert(0,'..')
    from memspectrum import MESA
except:
    from memspectrum import MESA

import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

if __name__ == "__main__":
    parser = OptionParser(usage)
    parser.add_option('--data', default=None, type='string', help='input data file')
    parser.add_option('--srate', default=None, type='float', help='sampling rate')
    parser.add_option('-T', default=None, type='float', help='duration of the data')
    (opts,args)=parser.parse_args()
    import time
    import matplotlib.pyplot as plt
    from scipy.signal import decimate
    
    srate = opts.srate
    dt = 1./srate
    T  = opts.T
    datafile = opts.data
    data = np.loadtxt(datafile)
    
    N = data.shape[0]
    f = np.fft.fftfreq(N, d=dt)
    t = np.arange(0,T,step=dt)
    M = MESA(data)
    start = time.perf_counter()
    P, ak, _ = M.solve(method = "Fast", optimisation_method = "FPE", m = int(2*N/(2*np.log(N))))
    print("p = {0}".format(len(ak)))
    elapsed = time.perf_counter()
    elapsed = elapsed - start
    print ("Time spent MESA: {0} s".format(elapsed))
    start = time.perf_counter()
    PSD    = M.spectrum(dt,f)
    elapsed = time.perf_counter()
    elapsed = elapsed - start
    print ("Time spent PSD: {0} s".format(elapsed))

    fig = plt.figure(1)
    ax  = fig.add_subplot(111)
    ax.loglog(f[:N//2], M.spectrum(dt,f)[:N//2],'-k')
    ax.set_xlim(1,srate/2.)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("PSD (Hz$^{-1}$)")
    plt.show()

