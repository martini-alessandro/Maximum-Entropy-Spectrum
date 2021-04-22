import sys
sys.path.insert(0,'..')
from memspectrum import MESA
from welch import *
import numpy as np
import matplotlib.pyplot as plt
#from init_plotting import init_plotting
from style_sheet import *
import pandas as pd
from scipy.interpolate import interp1d

plot = True
plot_same = True
compute = False
generate_fake_data = False
use_fake_data = True

plot_welch = True #if plot_same is true, this decides whether to plot welch method or MESA

	#folder to save the plot at
save_folder = '../paper/Images/comparison_LVC_data/'
true_PSD_file = 'GWTC1_GW150914_PSDs.dat'
#save_folder = 'comparison_LVC_data' #if on the cluster
srate = 4096.

T_list = [1, 5, 10,100, 1000] #list of times to make the comparison at
seglen = [512, 1024, 2048, 8192, 32768]

if generate_fake_data:
	PSD = np.loadtxt(true_PSD_file)
	freqs, PSD = PSD[:,0], PSD[:,1]
	import GenerateTimeSeries
	times, time_series, frequencies, frequency_series, psd_int = GenerateTimeSeries.generate_data(freqs, PSD, 1002., srate)#, np.min(freqs), np.max(freqs))
	np.savetxt('fake_data_4KHZ-1000.txt', time_series)
	
if compute:
	if use_fake_data:
		print("Using fake data")
		data = pd.read_csv("fake_data_4KHZ-1000.txt").to_numpy() #for fake data
	else:
		print("Using real data")
		#data = np.loadtxt("../examples/data/V-V1_GWOSC_4KHZ_R1-1186741846-32.txt")
		#data = pd.read_csv("../../GWAnomalyDetection/maxent/H-H1_GWOSC_16KHZ_R1-1126259447-32.txt.gz", skiprows = 3).to_numpy()
		data = pd.read_csv("H-H1_GWOSC_4KHZ_R1-1246525177-4096.txt.gz", skiprows = 3).to_numpy()

	data = np.squeeze(data)
	print("Loaded data: shape {}; srate {}; length {}s".format(data.shape, srate, len(data)/srate))

	for i, T in enumerate(T_list):
		print("Batch length: {}s".format(T))
		data_T = data[-int(srate*T)-4096*5:-4096*5]
		
		M = MESA()
		M.solve(data_T, early_stop = True, method = 'Standard')
		print("\tDone MESA")
		freqs, PSD_Welch = psd(data_T, srate, seglen[i]/float(srate),
			window_function  = None,
			overlap_fraction = 0.5,
			nfft = None,
			return_onesided = False)
		print("\tDone Welch")
		PSD_MESA = M.spectrum(1./srate, freqs)
		
		np.savetxt("plot_data/plot_{}_{}.txt".format(T, use_fake_data), np.column_stack([freqs, PSD_MESA, PSD_Welch]))

if plot_same:
	if use_fake_data:
		true_PSD = np.loadtxt(true_PSD_file)

	fig = init_plotting()
	fig.clf()
	fig.set_size_inches(3.4,3.4*2)
	ax = fig.gca()

	ax.set_xlabel(r"$f(Hz)$")
	ax.set_ylabel(r"$PSD (a.u.)$") #No ylabel
	
	yticks = []
	
	for i, T in enumerate(T_list):
		PSDs = np.loadtxt("plot_data/plot_{}_{}.txt".format(T, use_fake_data))
		N = PSDs.shape[0]
		freqs, PSD_MESA, PSD_Welch = PSDs[:int(N/2),0], PSDs[:int(N/2),1], PSDs[:int(N/2),2]

		offset = 1.
		if use_fake_data: step = 4.2
		else: step = 4.
		
			#scaling PSDs
		PSD_MESA =PSD_MESA/ 1e-41 *10**(step*i)*offset
		PSD_Welch = PSD_Welch/ 1e-41 * 10**(step*i)

		if plot_welch:
			filename = save_folder+"comparison_LVC_data_Welch.pdf".format(use_fake_data)
			ax.loglog(freqs, PSD_Welch , c = 'g', zorder = 0, label = "Welch", lw = 0.5)
		else:
			filename = save_folder+"comparison_LVC_data_MESA.pdf".format(use_fake_data)
			ax.loglog(freqs, PSD_MESA , c = 'r', zorder = 1, label = "MESA", lw = 0.5)
		
		yticks.append(10**(step*i)/1000.)
		print(yticks)
		
		if i == 0:
			ax.legend(loc = 'upper center', ncol = 1, fontsize = 10)
		
		if use_fake_data:
			ax.loglog(true_PSD[:,0], true_PSD[:,1]/1e-41*10**(step*i), '--', c = 'k', zorder = 2, lw = 0.75)
			#ax.loglog(true_PSD[:,0], true_PSD[:,1]/1e-41*10**(step*i)*offset, '--', c = 'k', zorder = 2, label = "True", lw = 0.75)
		
	ax.set_xlim(20,1024)
	ax.set_ylim(1e-8, 10**(step*len(T_list))/1000 )

	tick_list = ['T = {} s'.format(t) for t in T_list]
	ax.set_yticks(yticks)
	ax.set_yticklabels(tick_list, rotation = 'vertical', fontdict = {'verticalalignment':'center'})
	ax.yaxis.set_tick_params(width=0)

	plt.savefig(filename)
	print("Save file @ {}".format(filename))

	plt.show()
	quit()
	
if plot:
	if use_fake_data:
		true_PSD = np.loadtxt(true_PSD_file)
		
	fig_overall = init_plotting()
	fig_overall.clf()
	fig_overall.set_size_inches(3.4,3.4*2)
	for i, T in enumerate(T_list):
		PSDs = np.loadtxt("plot_data/plot_{}_{}.txt".format(T, use_fake_data))
		N = PSDs.shape[0]
		freqs, PSD_MESA, PSD_Welch = PSDs[:int(N/2),0], PSDs[:int(N/2),1], PSDs[:int(N/2),2]
		
		fig = init_plotting()
		ax = fig.gca()
		ax.set_title(r"$T = {}s$".format(T))
		ax.loglog(freqs, PSD_Welch, c = 'b', zorder = 0)
		ax.loglog(freqs, PSD_MESA, c = 'r', zorder = 1)
		if use_fake_data:
			ax.loglog(true_PSD[:,0], true_PSD[:,1], '--', c = 'k', zorder = 2)
		#ax.set_xlim(10,np.max(true_PSD[:,0]))
		ax.set_xlim(10,1024)
		ax.set_xlabel(r"$f(Hz)$")
		ax.set_ylabel(r"$PSD \left(\frac{1}{Hz} \right)$")
		#fig.subplots_adjust(left = 0.25, bottom = 0.25)
		plt.tight_layout(pad = 0.5)
		filename = save_folder+"comparison_LVC_data_T{}_fake_{}.pdf".format(T, use_fake_data)
		plt.savefig(filename)
		print("Save file @ {}".format(filename))
		
			#overall plot
		ax = fig_overall.add_subplot(5,1,i+1)#, sharex = True)
#		set_ax_style(ax)
		ax.set_title(r"$T = {}s$".format(T))
		ax.loglog(freqs, PSD_Welch, c = 'b', zorder = 0, label = "Welch")
		ax.loglog(freqs, PSD_MESA, c = 'r', zorder = 1, label = "MESA")
		if use_fake_data:
			ax.loglog(true_PSD[:,0], true_PSD[:,1], '--', c = 'k', zorder = 2, label = "True")
		#ax.set_ylabel(r"$PSD \left(\frac{1}{Hz} \right)$") #No ylabel
		ax.set_xlim(20,1024)
		if i == len(T_list)-1:
			ax.set_xlabel(r"$f(Hz)$")
		if i == 0:
			ax.legend(loc = 'upper center', ncol = 2+int(use_fake_data))#, fontsize='x-small')
	fig_overall.tight_layout()
	filename = save_folder+"comparison_LVC_data_overall_fake_{}.pdf".format(use_fake_data)
	plt.savefig(filename)
	print("Save file @ {}".format(filename))
	fig_overall.savefig(filename)
	
	plt.show()

		



















	
