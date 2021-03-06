#####
#Loading climate data
#Request ID: 72b941a9-4c27-4d17-8061-3525cfcbd646

#Variable:2m
#temperatureYear:2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020
#Month:January, February, March, April, May, June, July, August, September, October, November, December
#Day:01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
#Time:00:00, 01:00, 02:00, 03:00, 04:00, 05:00, 06:00, 07:00, 08:00, 09:00, 10:00, 11:00, 12:00, 13:00, 14:00, 15:00, 16:00, 17:00, 18:00, 19:00, 20:00, 21:00, 22:00, 23:00
#Sub-region extraction:North 45.5°, West 9.1°, South 45.4°, East 9.2°
#Format:GRIB

######

#PROBLEM: the code does not catch the trend over the year timescale!!

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'..')
import memspectrum
import scipy.signal

sys.path.insert(0,'../scripts')
import welch


#data = np.loadtxt('data_@1h_len95689.dat')[(195-30)*24:(195+30)*24] #summer
#data = np.loadtxt('data_@1h_len95689.dat')[(2567-30)*24:(2567+30)*24] #winter
#data = np.loadtxt('data_@1h_len95689.dat')
data = np.loadtxt('data_@1h_len271009_1989-2020.dat')

srate = 1./3600 # 1/(1 hours)
T = len(data)/srate

	#togliere ciclo annuale
	#togliere trend delle medie annuali

mu = np.mean(data)
f_0 = 1/(365)
t = np.linspace(0, len(data)/24., len(data))
A = (max(data)-min(data))/2.*0.5
phi_0 = 2.8
if True:
	def f(t, A, f_0,phi_0):
		return A*np.cos(2*np.pi*f_0*t +phi_0) +mu
	print("Initial guess for amplitude, freq and phase for annual trend: ",A, f_0,phi_0)
	popt, pcov = scipy.optimize.curve_fit(f, t, data, (A, f_0,phi_0))
	print("Final guess for amplitude, freq and phase for annual trend: ",*popt)
	A, f_0,phi_0 = popt
	y = f(t, A, f_0,phi_0)
else:
	y = f(t, A, f_0,phi_0)
#data = data -y

#computing spectrum
M = memspectrum.MESA()
if False:
	M.solve(data, method ='fast', optimisation_method = 'CAT', early_stop = False)
	M.save("model_CAT_long")
else:
	M.load("model_CAT_long")
N = len(data)
spec, f = M.spectrum(1./srate)
spec = spec[:int(N/2)]
f = f[:int(N/2)]

if False:
	f_welch, spec_welch = welch.psd(data, srate, len(data))
	spec = spec[:int(len(spec)/2)-1]
	f = f[:int(len(f)/2)-1]

#forecasting
N_forecast = 24*200 #1000 day of forecasting
predictions = M.forecast(data[:-N_forecast], N_forecast, 100)

#computing correlation
correlation = scipy.signal.correlate(data-np.mean(data), data-np.mean(data), mode="full")
#correlation = scipy.signal.correlate(data, data, mode="full")
correlation = correlation/max(correlation)
lags = scipy.signal.correlation_lags(data.size, data.size, mode="full")
	#doing fft
#spec_welch = np.fft.fft(correlation)
#plt.loglog(np.abs(spec_welch))


	#plot time series
fig = init_plotting()
ax = fig.gca()
t = np.linspace(0, len(data)/24., len(data))
plt.plot(t, data)
plt.plot(t, y)
plt.xlabel("Days")
plt.ylabel("Temperature (K)")
plt.savefig("temp_timeseries.pdf")

	#plot spectrum
plt.figure()
unit_shift = 3600*24 #frequency is made 1/day
plt.loglog(f * unit_shift , spec, label = "spectrum", c = 'k')
#plt.loglog(f_welch * unit_shift , spec_welch, label = "spectrum - welch", c = 'g')
plt.axvline(1./(3600.)*unit_shift, label = "hour", c = 'r', ls = '--')
plt.axvline(1./(3600.*24.)*unit_shift, label = "day", c = 'b', ls = '--')
plt.axvline(1./(3600.*12.)*unit_shift, label = "1/2 day", c = 'b', ls = '--')
plt.axvline(1./(3600.*24.*365.)*unit_shift, label = "year", c = 'g', ls = '--')
plt.xlabel("frequency (1/day)")
plt.ylabel("PSD (1/Hz)")
plt.legend()
plt.savefig("spectrum.pdf")

#plt.show()
#quit()

	#plot forecasting
plt.figure()
plt.title("{} days of forecasting".format(int(N_forecast/24)))
l, m, h = np.percentile(predictions,[5,50, 95],axis=0)
plt.plot(t[-(M.get_p()+N_forecast):],data[-(M.get_p()+N_forecast):],linewidth=1.5,color='r',zorder=2, label = 'data')
plt.axvline(t[-N_forecast])
plt.plot(t[-N_forecast:],m,'--',linewidth=0.5, color='k', zorder = 0)
plt.fill_between(t[-N_forecast:],l,h,facecolor='turquoise',alpha=0.5, zorder = 1)
plt.savefig("forecast.pdf")

	#plot a_k
plt.figure(figsize = (20,20))
N_start = int(len(correlation)/2)
delta_T = len(M.a_k)

correlation_spec = np.fft.irfft(spec*np.sqrt(np.mean(data))) #or there is a +1 in there...
correlation_spec/=correlation_spec[0]
#correlation_spec += np.square(np.mean(data))
#plt.plot(correlation_spec)
#plt.show()

plt.plot(lags[N_start:N_start+delta_T], correlation[N_start:N_start+delta_T], label='correlation')
plt.plot(lags[N_start:N_start+delta_T], correlation_spec[:delta_T], label = 'mesa correlation')
plt.legend()
np.savetxt('autocorr.dat', correlation[N_start:N_start+delta_T*10])

plt.show()

plt.plot(range(len(M.a_k)), M.a_k, 'o', ms =3, label = r'$a_k$')
for i in range(0,int(len(M.a_k)/24)):
	plt.axvline(i*24.+1, ls ='-.', lw = 1, c = 'g')
	plt.axvline(i*24., ls ='-.', lw = 1, c = 'r')
plt.xlabel('k')
plt.ylabel(r'$a_k$/autocorr')

N_start = int(len(correlation)/2)
delta_T = len(M.a_k)
plt.plot(lags[N_start:N_start+delta_T], (correlation[N_start:N_start+delta_T]-correlation[N_start+delta_T])/(correlation[N_start]-correlation[N_start+delta_T]), c = 'k', label = 'scaled_autocorr')

plt.legend()
plt.savefig("a_k.pdf")

plt.show()




