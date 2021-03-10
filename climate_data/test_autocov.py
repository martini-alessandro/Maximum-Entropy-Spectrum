"""
Two misteries:
	1) Why the empirical autocovariance and the mesa autocovariance have different scales?
	2) Which is the relation between autocovariance and autocorrelation?
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'..')
import memspectrum
import scipy.signal

sys.path.insert(0,'../scripts')
import welch

def empirical_autocorrelation(data):
	corr = scipy.signal.correlate(data, data, mode="full")
		#computing average value
	lags = scipy.signal.correlation_lags(data.size, data.size, mode="full")
	N_divide = len(data)-np.abs(lags) #number of data summed at each lag
	return np.divide(corr,N_divide)


data = np.loadtxt('data_@1h_len95689.dat')
#data = data - np.mean(data)# + 3.242

srate = 1/3600.

M = memspectrum.MESA()
M.solve(data, method ='fast', optimisation_method = 'FPE', early_stop = True)
M.save("test_model_FPE")
M.load("test_model_FPE")

K = -2.342932 #1- np.mean(data)*np.sum(M.a_k)
M.a_k = M.a_k/K
spec_scaled, f = M.spectrum(1./srate)

M.solve(data, method ='fast', optimisation_method = 'FPE', early_stop = True)
spec, f = M.spectrum(1./srate)
S_0 = spec[np.where(f==0)]
print(S_0)

	#spectrum of zero mean data
M_zero = memspectrum.MESA()
M_zero.solve(data- np.mean(data), method ='fast', optimisation_method = 'FPE', early_stop = True)
spec_zero, f = M_zero.spectrum(1./srate)
S_0_zero = spec_zero[np.where(f==0)]
print(S_0_zero)

print(np.allclose(spec_zero, spec), np.allclose(spec, spec_scaled*1/K**2))

plt.figure()
plt.loglog(f,spec_scaled*1/K**2, label = "symmetry")
plt.loglog(f,spec_zero, label = 'zero mean')
plt.loglog(f,spec, label = 'normal data')
plt.legend()
#plt.show()


autocov_empirical = empirical_autocorrelation(data-np.mean(data))
autocov = scipy.signal.correlate(data-np.mean(data), data-np.mean(data), mode="full")/len(data) #REMEMBER TO DIVIDE EVERYTHING BY N (it's an average value)
autocov /= max(autocov)
autocorr_empirical = empirical_autocorrelation(data)
autocorr_empirical -= np.square(np.mean(data))
autocorr = scipy.signal.correlate(data, data, mode="full")/len(data)
autocorr *= (np.var(data)+np.square(np.mean(data)))/max(autocorr)
autocorr -= np.square(np.mean(data))

lags = scipy.signal.correlation_lags(data.size, data.size, mode="full")

autocov_mesa_zero = M_zero.compute_autocovariance(3600., False)
autocov_mesa_zero *= np.var(data)/max(autocov_mesa_zero)
autocov_mesa = M.compute_autocovariance(3600., False)
autocov_mesa *= (np.var(data)+np.square(np.mean(data)))/max(autocov_mesa)
autocov_mesa -= np.square(np.mean(data))

N_start = int(len(autocov)/2)
delta_T = 5000

print("Autocorr values",autocorr_empirical[N_start:N_start+5], autocov_empirical[N_start:N_start+5])

print(autocorr[N_start], S_0, autocov[N_start], S_0_zero)
print(autocorr[N_start]/autocov[N_start], (np.mean(data)))
print(autocorr/autocov)
plt.figure()
plt.plot(autocorr-np.mean(data)**2)
plt.plot(autocov)

plt.figure()
#plt.plot(lags[N_start:N_start+delta_T], autocov_mesa_zero[:delta_T], label = "MESA zero")
plt.plot(lags[N_start:N_start+delta_T], autocov_mesa[:delta_T], label = "MESA")
plt.plot(lags[N_start:N_start+delta_T], autocorr[N_start:N_start+delta_T], label = "autocov")
#plt.plot(lags[N_start:N_start+delta_T], autocov_empirical[N_start:N_start+delta_T], label = "autocov empirical")
#plt.plot(lags[N_start:N_start+delta_T], autocorr_empirical[N_start:N_start+delta_T], label = "autocorr empirical")
#plt.yscale('log')
plt.legend()
plt.show()

