"""
Script to generate LIGO noise similar to the design psd at a fixed autoregressive order p.

To generate 32 s of data @ 4096 Hz with AR order 300 run

131073


"""
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
import memspectrum
from memspectrum.GenerateTimeSeries import generate_data

import argparse
import warnings


#############################################

parser = argparse.ArgumentParser(__doc__)
parser.add_argument(
	"--p", default = 200, type = int, 
	help="Autoregressive order of the final timeseries")
parser.add_argument(
	"--srate", default = 4096, type = float, 
	help="Sample rate of the final timeseries")
parser.add_argument(
	"--t", default = 32, type = float, 
	help="Length in seconds of the final timeseries")
parser.add_argument(
	"--savefile", default = None, type = str, 
	help="File to save the timeseries to")

args, filenames = parser.parse_known_args()

#############################################

psd_file = 'aligo_O3actual_H1.txt'
if not os.path.isfile(psd_file):
	subprocess.run('wget https://dcc.ligo.org/public/0165/T2000012/002/aligo_O3actual_H1.txt', shell = True)

f_input, psd_input = np.loadtxt(psd_file).T

t_grid, time_series, f_grid, frequency_series, interpolated_psd = generate_data(f_input, psd_input, args.t, args.srate, seed = 0)

m = memspectrum.MESA()
m.solve(time_series, method = 'Standard', optimisation_method = 'Fixed', m = args.p+1)
assert m.p == args.p

f_mesa, psd_mesa = m.spectrum(dt = 1./args.srate, onesided = True)

burn_in = 20_000

fake_data = m.forecast(np.random.normal(0, m.P, args.p), length = int(args.t*args.srate)+burn_in, number_of_simulations = 1, P = None, include_data = False, seed = None, verbose = True)[0,burn_in:]

	#Checking if generate data are sane
m_check = memspectrum.MESA()
m_check.solve(fake_data, method = 'Standard')#, optimisation_method = 'Fixed', m = args.p+1)
f_check, psd_check = m_check.spectrum(dt = 1./args.srate, onesided = True)
print("True data p: ", m.p)
print("Generated data p: ", m_check.p)
if np.abs(m_check.p - m.p)>1:
	warnings.warn("The PSD of the fake data is not super sane...")

if args.savefile:
	head = 'T = {} | srate = {} | p = {}'.format(args.t, args.srate, args.p)
	np.savetxt(args.savefile, fake_data, header = head)

	#Plotting

plt.figure()
plt.loglog(f_input, psd_input, ls = '--', label = 'true PSD')
plt.loglog(f_mesa, psd_mesa, label = 'estimated PSD')
plt.loglog(f_check, psd_check, label = 'fake data PSD')
plt.xlabel('f (Hz)')
plt.ylabel('PSD (1/Hz)')
plt.legend()
plt.title('Estimated PSD (vs true PSD)')

plt.show()
















