import numpy as np

try:
	import sys
	sys.path.insert(0,'..')
	from memspectrum import MESA
	from GenerateTimeSeries import generate_data
except:
	from memspectrum import MESA
	from memspectrum.GenerateTimeSeries import generate_data


import matplotlib.pyplot as plt

def lisa_psd(f):
    # see https://arxiv.org/pdf/1201.3621.pdf
    L       = 1e9
    C       = 3e8
    psd     = np.zeros(f.shape[0])
    Sacc    = (1.37e-32)*(1+1e-4/f)/f**4
    Sxsn    = 5.25e-23
    Sxomn   = 6.28e-23
    psd     = (20./3.)*(4*Sacc+Sxsn+Sxomn)/L**2
    psd    *= (1.0+(f/(0.41*0.5*C/L))**2)
    
    return psd
if __name__ == "__main__":
    srate = 1.0
    dt = 1./srate
    bandpassing = 0
    fmin = 1e-5
    fmax = srate/2.0
    T  = 1e5
    df = 1./T
    f  = np.arange(fmin, fmax, step = df)
    psd = lisa_psd(f)
    t, time_series, f, frequency_series, psd = generate_data(f,
                          psd,
                          T = T,
                          sampling_rate = srate,
                          fmin = fmin,
                          fmax = fmax,
                          zero_noise = False,
                          asd = False)
    N = time_series.shape[0]
    f = np.fft.fftfreq(N, d=dt)
    M = MESA(time_series)
    P, ak, _ = M.solve(method = "Fast", optimisation_method = "FPE", m = int(2*N/(2*np.log(N))))
    print("p = {0}".format(len(ak)))
    PSD     = M.spectrum(dt,f)

    fig = plt.figure(1)
    ax  = fig.add_subplot(111)
    ax.loglog(f[:N//2], M.spectrum(dt,f)[:N//2],'-k')
    ax.loglog(f[:N//2], lisa_psd(f[:N//2]),'--r')
    ax.set_xlim(fmin,fmax)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("PSD (Hz$^{-1}$)")
    
    M = MESA(time_series[:int(0.9*N)])
    P, ak, _ = M.solve(method = "Fast", optimisation_method = "FPE", m = int(2*N/(2*np.log(N))))
    Np = 3
    prediction = M.forecast(int(0.1*N), Np)
    l, h = np.percentile(prediction,[5,95],axis=0)
    fig = plt.figure(2)
    ax  = fig.add_subplot(111)
    ax.plot(t,time_series,linewidth=1.5,color='r',zorder=2)
    ax.axvline(t[int(0.9*N)],linestyle='dashed',color='blue')
    for p in range(Np):
        ax.plot(t[int(0.9*N):],prediction[p,:],linewidth=0.5, color='k', zorder = 0)
    ax.fill_between(t[int(0.9*N):],l,h,facecolor='turquoise',alpha=0.5, zorder = 1)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("strain")
    plt.show()

