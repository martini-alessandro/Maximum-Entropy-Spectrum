"""
This script load a audio file holding noise and it computes its PSD. It also produces another audio file with syntetic noise, sharing the same features as the original noise file.
"""

import numpy as np
import scipy.io.wavfile
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'..')
from mesa.mesa import *
import mesa.GenerateTimeSeries

t = 4. #seconds of data

	#loading data and preparing input to MESA
rate, data = scipy.io.wavfile.read("waterfall_data/waterfall_noise.wav") #loading waterfall noise file
	#data is (N,2): stereophonic sound
data_MESA = data[:int(t*rate),0].astype(np.float64)
#data_MESA = data_MESA+1j*0. #why complex data?
dt = 1./rate
times_out = np.linspace(0., len(data_MESA)*dt, len(data_MESA))

	#computing PSD with MESA
M = MESA(data_MESA)
P, ak, opt= M.solve(method = "Fast", optimisation_method = "FPE", m = int(2*len(data_MESA)/(2*np.log(len(data_MESA)))))

	#evaluating the spectrum
N_points = 100000
f = np.linspace(0,10000,N_points) #Why this?? Is only employed for getting the N mesa.spectrum 
PSD =  M.spectrum(dt,len(f))[0][:int(N_points/2)] #we want only positive frequencies...
f_PSD = np.fft.fftfreq(len(f),dt)[:int(N_points/2)] #this it the actual frequency grid that PSD is evaluated at

	#evaluating the spectrum bis
f_bis = np.linspace(0,20000,10000)
PSD_bis =  M.spectrum_bis(f_bis, dt)

	#generate syntetic noise
times, time_series, frequencies, frequency_series, psd_int = mesa.GenerateTimeSeries.generate_data(f_PSD, PSD.real, T=16., sampling_rate = rate, zero_noise = False)

	#saving noise to file
to_save = np.column_stack([time_series, time_series])
scipy.io.wavfile.write("waterfall_data/simulated_noise.wav", rate, to_save.real.astype(np.int16))

	#checking whether psd is correct in the simulated noise
	#simulated noise is loaded and psd is estimated again
rate, data = scipy.io.wavfile.read("waterfall_data/simulated_noise.wav") #loading waterfall noise file
data_MESA = data[:int(t*rate),0].astype(np.float64)
data_MESA = data_MESA+1j*0.
dt = 1./rate
M = MESA(data_MESA)
P, ak, opt= M.solve(method = "Fast", optimisation_method = "FPE", m = int(2*len(data_MESA)/(2*np.log(len(data_MESA)))))
PSD_sim =  M.spectrum(dt,len(f))[0][:int(N_points/2)]

print("All done: if you like to listen to the output file, type \"aplay waterfall_data/simulated_noise.wav\" ")
#plotting some quantities

#frequency series of the reconstructed data
fig, ax = plt.subplots(2,1, sharex = True)
plt.suptitle("Real part of the frequency series of the simulated data with the PSD")
ax[0].plot(frequencies,frequency_series.real)
ax[0].plot(frequencies, np.sqrt(psd_int), label = "sqrt(PSD)")
ax[0].legend()

ax[1].plot(frequencies, psd_int)
ax[1].set_yscale('log')
ax[1].set_xlabel("frequency (Hz)")


#comparison between original psd and "simulated" psd
plt.figure()
plt.title("Comparison between PSD computed from original data (empirical)\nand simulated data (simulated and bis)")
plt.plot(f_PSD, PSD.real, label = "empirical")
plt.plot(f_PSD, PSD_sim.real, label = "simulated")
plt.plot(f_bis, PSD_bis, label = "bis", ms = 100)
plt.yscale('log')
plt.xlabel("frequency (Hz)")
plt.legend()

#comparison between two methods for spectrum
plt.figure()
plt.title("Difference between spectrum and spectrum bis")
PSD_std = np.interp(f_bis, f_PSD, PSD.real)
plt.plot(f_bis, (PSD_bis - PSD_std)/PSD_std)
plt.xlabel("frequency (Hz)")
plt.show()





