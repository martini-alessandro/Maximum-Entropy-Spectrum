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

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'..')
import memspectrum

data = np.loadtxt('data_@1h_len95689.dat')

srate = 1./3600 # 1/(1 hours)
T = len(data)/srate

#computing spectrum
M = memspectrum.MESA()
M.solve(data)
f = np.linspace(1./(T), 0.5*srate,1000)
spec = M.spectrum(1./srate, f)

plt.figure()
plt.plot(np.linspace(0, len(data)/24., len(data)), data)
plt.xlabel("Days")
plt.ylabel("Temperature (K)")
plt.savefig("temp_timeseries.pdf")

plt.figure()
unit_shift = 3600*24 #frequency is made 1/day
plt.loglog(f * unit_shift , spec, label = "spectrum", c = 'k')
plt.axvline(1./(3600.)*unit_shift, label = "hour", c = 'r', ls = '--')
plt.axvline(1./(3600.*24.)*unit_shift, label = "day", c = 'b', ls = '--')
plt.axvline(1./(3600.*24.*365.)*unit_shift, label = "year", c = 'g', ls = '--')
plt.xlabel("frequency (1/day)")
plt.ylabel("PSD")
plt.legend()
plt.savefig("spectrum.pdf")

plt.show()
