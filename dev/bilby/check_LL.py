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

GPS_time_list = np.array([1238170000+i*1000 for i in range(100)])
duration = 16
sampling_frequency = 4096
ifo_str = 'L1'
channel_name = '{}:GDS-CALIB_STRAIN'.format(ifo_str)
minimum_frequency, maximum_frequency = 20, 1024
load = True

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
if not load:
	ll_list = []
	for GPS_time in GPS_time_list:
		#creating the ifos and loading the data with gwpy (function set_strain_data_from_channel_name)
		ifos = bilby.gw.detector.InterferometerList([ifo_str])
		for ifo in ifos:
			ifo.minimum_frequency = minimum_frequency
			ifo.maximum_frequency = maximum_frequency
				#This is a dirty trick to make package nds2 (required by gwpy) available: works only on CIT
				#How to install this package in a normal way? I haven't figured out yet :(
			#sys.path.append('/usr/lib64/python3.6/site-packages')
			#sys.path.append('/usr/lib/python3.6/site-packages')
			ifo.set_strain_data_from_channel_name(channel_name,
				sampling_frequency = sampling_frequency, duration = duration, start_time= GPS_time)
			#sys.path.remove('/usr/lib64/python3.6/site-packages')
			#sys.path.remove('/usr/lib/python3.6/site-packages')

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
			mesa_ll = np.nan
			std_ll = np.nan

		print("GPS time: ", GPS_time)
		print("\tMesa ll: ", mesa_ll)
		print("\tStandard ll: ", std_ll)
		
		ll_list.append((mesa_ll, std_ll))

	ll_list = np.array(ll_list) #(N_GPS, 2 )
	np.savetxt('plots/ll_list.dat', ll_list)
else:
	ll_list = np.loadtxt('plots/ll_list.dat')

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

plt.figure()
plt.title('Cumulative LL diff')
LL_cum_diff = ll_list[:,0]-ll_list[:,1]
ids_nan = np.isnan(LL_cum_diff)
LL_cum_diff = np.nancumsum(LL_cum_diff)
LL_cum_diff[ids_nan] = np.nan

plt.plot(GPS_time_list - GPS_time_list[0], LL_cum_diff, 'o', label = 'duration = {}s'.format(duration))
plt.legend()
plt.xlabel("Time since GPS {} (s)".format(GPS_time_list[0]))
plt.savefig("plots/LL_cum_diff_plot.png")

plt.figure()
LL_diff = (ll_list[:,0]-ll_list[:,1])#/1e6
LL_diff = LL_diff[~ids_nan]
plt.title('LL diff cumualtive hist')
plt.hist(np.log10(LL_diff), bins = 300, cumulative = True, histtype='step')
plt.xlabel(r"$\log_{10}(LL_{MESA}-LL_{std})$")
plt.savefig("plots/LL_diff_cum_hist.png")

plt.figure()
percentile = 90
plt.title('LL diff hist (<{}%)'.format(percentile))
LL_diff = LL_diff[np.where(LL_diff<np.percentile(LL_diff, percentile))]
plt.hist(np.log10(LL_diff), bins = 10, histtype='step')
plt.xlabel(r"$\log_{10}(LL_{MESA}-LL_{std})$")
plt.savefig("plots/LL_diff_hist.png")

plt.show()


















