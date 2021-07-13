try:
	import sys
	sys.path.insert(0,'..')
	from memspectrum import MESA
except:
	from memspectrum import MESA

import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    import time
    import matplotlib.pyplot as plt
    from scipy.signal import decimate
    
    srate = 2048.
    dt = 1./srate
    bandpassing = 0
    f_min_bp = 20.0
    f_max_bp = (srate-20)/2.0
    T  = 10
    datafile = "data/V-V1_GWOSC_4KHZ_R1-1186741846-32.txt"
    data = np.loadtxt(datafile)[:int(T*4096)]
    if srate != 4096.:
        data = decimate(data, int(4096/srate), zero_phase=True)
    if bandpassing == 1:
        from scipy.signal      import butter, filtfilt
        bb, ab = butter(4, [f_min_bp/(0.5*srate), f_max_bp/(0.5*srate)], btype='band')
        data = filtfilt(bb, ab, data)
    
    N = data.shape[0]
    f = np.fft.fftfreq(N, d=dt)
    t = np.arange(0,T,step=dt)
    M = MESA()
    start = time.perf_counter()
    P, ak, _ = M.solve(data,method = "Fast", optimisation_method = "FPE", m = int(2*N/(2*np.log(N))))
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
    ax.loglog(f[:N//2], M.spectrum(dt,f)[:N//2],'-k', label = 'm = {}'.format(len(ak)))
    M.solve(data,method = "Fast", optimisation_method = "FPE", m =100)
    ax.loglog(f[:N//2], M.spectrum(dt,f)[:N//2],'-r', label = 'm = 100')
    ax.set_xlim(1,srate/2.)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("PSD (Hz$^{-1}$)")
    ax.legend()
    
    M = MESA()
    P, ak, _ = M.solve(data[:int(0.75*N)],method = "Fast", optimisation_method = "FPE", m = int(2*N/(2*np.log(N))))
    Np = 100
    prediction = M.forecast(data[:int(0.75*N)], int(0.25*N), Np)
    l, h = np.percentile(prediction,[5,95],axis=0)
    fig = plt.figure(2)
    ax  = fig.add_subplot(111)
    ax.plot(t,data,linewidth=1.5,color='r',zorder=2)
    ax.axvline(t[int(0.75*N)],linestyle='dashed',color='blue')
    ax.plot(t[int(0.75*N):],prediction.T,linewidth=0.5, color='k', zorder = 0)
    ax.fill_between(t[int(0.75*N):],l,h,facecolor='turquoise',alpha=0.5, zorder = 1)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("strain")
    plt.show()
