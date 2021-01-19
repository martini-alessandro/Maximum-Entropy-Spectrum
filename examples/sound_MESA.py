"""
This script load a audio file holding noise and it computes its PSD. It also produces another audio file with syntetic noise, sharing the same features as the original noise file.
"""

import numpy as np
import scipy.io.wavfile
import matplotlib.pyplot as plt

try:
    import sys
    sys.path.insert(0,'..')
    from mesa import *
    import GenerateTimeSeries
except:
    from mesa import *
    import mesa.GenerateTimeSeries as GenerateTimeSeries

t = 4. #seconds of data

    #loading data and preparing input to MESA
rate, data = scipy.io.wavfile.read("data/waterfall_noise.wav") #loading waterfall noise file
    #data is (N,2): stereophonic sound
data_MESA = data[:int(t*rate),0].astype(np.float64)
dt = 1./rate
times_out = np.linspace(0., len(data_MESA)*dt, len(data_MESA))

    #computing PSD with MESA
M = MESA(data_MESA)
P, ak, opt= M.solve(method = "Fast", optimisation_method = "FPE", m = int(2*len(data_MESA)/(2*np.log(len(data_MESA)))))

    #evaluating the spectrum
N_points = 100000
f_PSD = np.linspace(0,20000,N_points) 
PSD = M.spectrum(dt, f_PSD)

    #generate syntetic noise
times, time_series, frequencies, frequency_series, psd_int = GenerateTimeSeries.generate_noise_mesa(M, T=16., sampling_rate = rate, N_series = 1)

    #saving noise to file
to_save = np.column_stack([time_series, time_series])
scipy.io.wavfile.write("data/simulated_noise.wav", rate, to_save.real.astype(np.int16))

    #checking whether psd is correct in the simulated noise
    #simulated noise is loaded and psd is estimated again
rate, data = scipy.io.wavfile.read("data/simulated_noise.wav") #loading waterfall noise file
data_MESA = data[:int(t*rate),0].astype(np.float64)
data_MESA = data_MESA+1j*0.
dt = 1./rate
M = MESA(data_MESA)
P, ak, opt= M.solve(method = "Fast", optimisation_method = "FPE", m = int(2*len(data_MESA)/(2*np.log(len(data_MESA)))))
PSD_sim =  M.spectrum(dt,f_PSD) #PSD computed on simulated data

print("All done: if you like to listen to the output file, type \"aplay data/simulated_noise.wav\" ")
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
plt.yscale('log')
plt.xlabel("frequency (Hz)")
plt.legend()

plt.show()





