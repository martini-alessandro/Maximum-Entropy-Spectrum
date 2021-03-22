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
	return np.divide(corr,N_divide) #*np.mean(np.square(data)
	

data = np.loadtxt('data_@1h_len95689.dat')
#data = data - np.mean(data)# + 3.242

srate = 1/3600.

M = memspectrum.MESA()
M.solve(data, method ='fast', optimisation_method = 'FPE', early_stop = True)
M.save("test_model_FPE")
M.load("test_model_FPE")

K = 2000.342932 #1- np.mean(data)*np.sum(M.a_k)
M.a_k = M.a_k/K
M.P = M.P/K**2
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
print(np.sum(M_zero.a_k), np.sum(M.a_k))

print(np.allclose(spec_zero, spec), np.allclose(spec, spec_scaled))

print(M.a_k, M_zero.a_k)
print(M.P, M_zero.P)
print(np.sum(M.a_k), np.sum(M_zero.a_k))
print(np.allclose(M.a_k[:10],M_zero.a_k[:10]))
plt.plot(M_zero.a_k)
plt.plot(M.a_k)
#plt.show()

#quit()

plt.figure()
plt.loglog(f,spec_scaled*1/K**0, label = "symmetry")
plt.loglog(f,spec_zero, label = 'zero mean')
plt.loglog(f,spec, label = 'normal data')
plt.legend()
plt.show()


autocov_empirical = empirical_autocorrelation(data-np.mean(data))
#autocov_empirical /= max(autocov_empirical)

autocorr_empirical = empirical_autocorrelation(data)
#autocorr_empirical -= np.square(np.mean(data))

if False:
	mu_avg_sq = np.zeros(autocorr_empirical.shape)
	mu = np.mean(data)
	for i in range(1,len(data)-1):
		mu_avg_sq[len(data)+i] = -np.mean(data[i:]+data[:-i]) * mu + mu**2
		mu_avg_sq[len(data)-i] = mu_avg_sq[len(data)+i]
	np.savetxt('mu_avg_sq.dat', mu_avg_sq)
else:
	mu_avg_sq = np.loadtxt('mu_avg_sq.dat')

mu_avg_sq[len(data)] = -np.square(np.mean(data))
#print("Mean: ", -np.mean(data[k:]+data[:-k]) * np.mean(data) +np.mean(data)**2, -np.mean(data)**2)

autocorr_empirical += mu_avg_sq

lags = scipy.signal.correlation_lags(data.size, data.size, mode="full")

autocov_mesa_zero = M_zero.compute_autocovariance(3600., False)
autocov_mesa_zero *= np.var(data)/max(autocov_mesa_zero)

autocov_mesa = M.compute_autocovariance(3600., True)
#autocov_mesa -= max(autocov_mesa)-np.mean(np.square(data))
#factor = max(autocov_mesa)/np.mean(np.square(data))
#autocov_mesa += mu_avg_sq[len(data):len(data)+len(autocov_mesa)]
#autocov_mesa *= max(autocov_mesa_zero)/max(autocov_mesa)
#print(autocov_mesa[0], max(autocov_mesa_zero))

N_start = len(data)
delta_T = 5000

plt.figure()
plt.plot(lags[N_start:N_start+delta_T], autocov_mesa[:delta_T])

plt.figure()
plt.plot(lags[N_start:N_start+delta_T], autocov_mesa_zero[:delta_T]/autocov_mesa_zero[0])
#plt.plot(lags[N_start:N_start+delta_T], autocorr_empirical[N_start:N_start+delta_T]-mu_avg_sq[N_start:N_start+delta_T])

plt.figure()
#plt.plot(lags[N_start:N_start+delta_T], autocov_mesa[:delta_T], label = "MESA")
plt.plot(lags[N_start:N_start+delta_T], autocov_mesa_zero[:delta_T], label = "MESA zero")
plt.plot(lags[N_start:N_start+delta_T], autocov_empirical[N_start:N_start+delta_T], label = "autocov empirical")
plt.plot(lags[N_start:N_start+delta_T], autocorr_empirical[N_start:N_start+delta_T], label = "autocorr empirical")
#plt.yscale('log')
plt.legend()
plt.show()

