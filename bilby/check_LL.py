from bilby_mesa import MESAGravitationalWaveTransient

import argparse

import numpy as np
import matplotlib.pyplot as plt
import bilby
from bilby.gw.likelihood import Likelihood #needs version 1.1.2 (the latest doesn't work, for some reason)
import memspectrum
from scipy.signal.windows import tukey
import sys

#####################

GPS_time_list = np.array([1238170000+i*1000 for i in range(50)])
duration = 16
sampling_frequency = 4096
ifo_str = 'L1'
channel_name = '{}:GDS-CALIB_STRAIN'.format(ifo_str)
minimum_frequency, maximum_frequency = 20, 1024

#mesa settings
N = duration*sampling_frequency
mesa_m = int(2*N/np.log(2.*N))
mesa_method = 'Fast'
mesa_optimisation_method = 'FPE'


#creating a WF generator (even though we won't need to generate WFs)
waveform_arguments = dict(
		waveform_approximant="TaylorF2",
		reference_frequency=minimum_frequency,
		minimum_frequency=minimum_frequency,
		catch_waveform_errors=True,
	)

waveform_generator = bilby.gw.WaveformGenerator(
		duration=duration,
		sampling_frequency=sampling_frequency,
		frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
		parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
		waveform_arguments=waveform_arguments,
	)



########
# Loops on several GPS times
########

ll_list = []
for GPS_time in GPS_time_list:
	#creating the ifos and loading the data with gwpy (function set_strain_data_from_channel_name)
	ifos = bilby.gw.detector.InterferometerList([ifo_str])
	for ifo in ifos:
		ifo.minimum_frequency = minimum_frequency
		ifo.maximum_frequency = maximum_frequency
			#This is a dirty trick to make package nds2 (required by gwpy) available: works only on CIT
			#How to install this package in a normal way? I haven't figured out yet :(
		sys.path.append('/usr/lib64/python3.6/site-packages')
		sys.path.append('/usr/lib/python3.6/site-packages')
		ifo.set_strain_data_from_channel_name(channel_name,
			sampling_frequency = sampling_frequency, duration = duration, start_time= GPS_time)
		sys.path.remove('/usr/lib64/python3.6/site-packages')
		sys.path.remove('/usr/lib/python3.6/site-packages')

		#plotting ifo data
	if np.any(np.abs(ifo.strain_data.time_domain_strain)>1e-17): bad = True
	else: bad = False
	
	plt.figure()
	if bad: plt.title("Bad strain!! Whyyyy?")
	plt.plot(np.linspace(0, duration, len(ifo.strain_data.time_domain_strain)), ifo.strain_data.time_domain_strain)
	plt.xlabel("Time since GPS {} (s)".format(GPS_time))
	plt.savefig('plots/strain_{}.png'.format(GPS_time))
	plt.close('all')
	

	standard_likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
				interferometers=ifos, waveform_generator=waveform_generator
			)

	mesa_likelihood = MESAGravitationalWaveTransient(
				interferometers=ifos,
				waveform_generator=waveform_generator,
				mesa_m=mesa_m,
				mesa_method=mesa_method,
				mesa_optimisation_method=mesa_optimisation_method,
			)
	mesa_likelihood.parameters["mesa_m"] = mesa_m

	try:
		assert not bad
		mesa_ll = mesa_likelihood.noise_log_likelihood()
		std_ll = standard_likelihood.noise_log_likelihood()
	except:
		print("Something wrong happened: bad data??")
		mesa_ll = 1.
		std_ll = 1.

	print("GPS time: ", GPS_time)
	print("\tMesa ll: ", mesa_ll)
	print("\tStandard ll: ", std_ll)
	
	ll_list.append((mesa_ll, std_ll))

ll_list = np.array(ll_list) #(N_GPS, 2 )


########
# Plotting
########

plt.figure()
plt.title('MESA LL/STD LL')
plt.plot(GPS_time_list - GPS_time_list[0], ll_list[:,0]/ll_list[:,1], 'o', label = 'duration = {}s'.format(duration))
plt.legend()
plt.xlabel("Time since GPS {} (s)".format(GPS_time_list[0]))
plt.savefig("plots/LL_ratio_plot.png")

plt.figure()
plt.title('MESA LL - STD LL')
plt.plot(GPS_time_list - GPS_time_list[0], ll_list[:,0]-ll_list[:,1], 'o', label = 'duration = {}s'.format(duration))
plt.legend()
plt.xlabel("Time since GPS {} (s)".format(GPS_time_list[0]))
plt.savefig("plots/LL_diff_plot.png")






















