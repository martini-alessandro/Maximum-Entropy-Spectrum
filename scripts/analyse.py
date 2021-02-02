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
    parser = OptionParser()
    parser.add_option('--data', default=None, type='string', help='input data file')
    parser.add_option('--srate', default=None, type='float', help='sampling rate')
    parser.add_option('-T', default=None, type='float', help='duration of the data')
    (opts,args)=parser.parse_args()
    import time
    import matplotlib.pyplot as plt
    from scipy.signal import decimate
    
    srate = opts.srate
    dt = 1./srate
    datafile = opts.data
    d = np.genfromtxt(datafile, delimiter=',', names=True)
    data = d['CloseLast'][::-1]
    T = len(data)
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
    print(PSD)
    
    elapsed = time.perf_counter()
    elapsed = elapsed - start
    print ("Time spent PSD: {0} s".format(elapsed))

    fig = plt.figure(1)
    ax  = fig.add_subplot(111)
    ax.loglog(f[:N//2], M.spectrum(dt,f)[:N//2],'-k')
#    ax.set_xlim(1,srate/2.)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("PSD (Hz$^{-1}$)")
    
    
    t = np.arange(T-1,T+199,step=dt)
    prediction = M.forecast(200, 1000)
    l, m, h = np.percentile(prediction,[5,50,95],axis=0)
    fig = plt.figure(2)
    ax  = fig.add_subplot(111)
    ax.plot(np.arange(0,T,step=dt),data,linewidth=1.5,color='r',zorder=2)
    ax.plot(t,prediction.T,linewidth=0.5, color='k', zorder = 0)
    ax.plot(t,m,linewidth=1.5, color='k', zorder=3)
    ax.fill_between(t,l,h,facecolor='turquoise',alpha=0.5, zorder = 1)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("strain")
    plt.show()
