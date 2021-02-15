import sys
sys.path.insert(0,'..')
from memspectrum import MESA
from welch import *
import numpy as np
import matplotlib.pyplot as plt
from init_plotting import init_plotting
import pandas as pd


plot = True
compute = False
generate_fake_data = False
use_fake_data = True

	#folder to save the plot at
save_folder = '../paper/Images/comparison_LVC_data/'
true_PSD_file = 'GWTC1_GW150914_PSDs.dat'
#save_folder = 'comparison_LVC_data' #if on the cluster
srate = 4096.

T_list = [1,5, 10,100, 1000] #list of times to make the comparison at

if generate_fake_data:
	PSD = np.loadtxt(true_PSD_file)
	freqs, PSD = PSD[:,0], PSD[:,1]
	import GenerateTimeSeries
	times, time_series, frequencies, frequency_series, psd_int =GenerateTimeSeries.generate_data(freqs, PSD, 1002., srate, np.min(freqs), np.max(freqs))
	np.savetxt('fake_data_4KHZ-1000.txt', time_series)
	
if compute:
	if use_fake_data:
		print("Using fake data")
		data = pd.read_csv("fake_data_4KHZ-1000.txt").to_numpy() #for fake data
	else:
		#data = np.loadtxt("../examples/data/V-V1_GWOSC_4KHZ_R1-1186741846-32.txt")
		#data = pd.read_csv("../../GWAnomalyDetection/maxent/H-H1_GWOSC_16KHZ_R1-1126259447-32.txt.gz", skiprows = 3).to_numpy()
		data = pd.read_csv("H-H1_GWOSC_4KHZ_R1-1247614487-4096.txt.gz", skiprows = 3).to_numpy()

	data = np.squeeze(data)
	print("Loaded data: shape {}; srate {}; length {}s".format(data.shape,srate, len(data)/srate))

	for i, T in enumerate(T_list):
		print("Batch length: {}s".format(T))
		data_T = data[:int(srate*T)]

		seglen = min(20.,T)
		freqs, PSD_Welch = psd(data_T, srate, seglen,
			window_function  = None,
			overlap_fraction = 0.5,
			nfft = None,
			return_onesided = True)
		M = MESA()
		M.solve(data_T, early_stop = True)
		PSD_MESA = M.spectrum(1./srate, freqs)
		
		np.savetxt("plot_data/plot_{}_{}.txt".format(T, use_fake_data), np.column_stack([freqs, PSD_MESA, PSD_Welch]))
		
		#plt.loglog(freqs, PSD_Welch)
		#plt.loglog(freqs, PSD_MESA)
		#plt.show()
	
	
if plot:
	if use_fake_data:
		true_PSD = np.loadtxt(true_PSD_file)
	for i, T in enumerate(T_list):
		PSDs = np.loadtxt("plot_data/plot_{}_{}.txt".format(T, use_fake_data))
		freqs, PSD_MESA, PSD_Welch = PSDs[:,0], PSDs[:,1], PSDs[:,2]
		
		fig = init_plotting()
		ax = fig.gca()
		ax.set_title(r"$T = {}s$".format(T))
		ax.loglog(freqs, PSD_Welch, c = 'b', zorder = 0)
		ax.loglog(freqs, PSD_MESA, c = 'r', zorder = 1)
		if use_fake_data:
			ax.loglog(true_PSD[:,0], true_PSD[:,1], '--', c = 'k', zorder = 2)
		ax.set_xlabel(r"$f(Hz)$")
		ax.set_ylabel(r"$PSD \left(\frac{1}{\sqrt{Hz}} \right)$")
		#fig.subplots_adjust(left = 0.25, bottom = 0.25)
		plt.tight_layout(pad = 0.5)
		plt.savefig(save_folder+"comparison_LVC_data_T{}_fake_{}.pdf".format(T, use_fake_data))
		plt.show()

		



















	
