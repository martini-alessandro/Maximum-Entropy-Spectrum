# https://github.com/ecmwf/cfgrib
# TO download data: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

ds = xr.open_dataset('adaptor.mars.internal-1614589133.4434555-10786-9-72b941a9-4c27-4d17-8061-3525cfcbd646.grib', engine='cfgrib', backend_kwargs={'read_keys': ['experimentVersionNumber']})
print(ds)
print(ds.t2m.attrs)
#help(xr.open_dataset)
#help(ds)

print("UFFAAAAAA")

print(ds.data_vars['t2m'])

print("UFFAAAAAAAAAAAAAAAAAAAAA\n")

print(ds.data_vars['t2m'].time) #actually keeps the day
print(ds.data_vars['t2m'].step)
print(ds.data_vars['t2m'].latitude)
print(ds.data_vars['t2m'].longitude)
print(ds.data_vars['t2m'].shape)

data = ds.data_vars['t2m']

	#doing a time series
b = np.array(ds.data_vars['t2m'][:,:,0,0])
print(b.shape)
print(b[:,:])

b = b.flatten()
print(b)
b = b[~np.isnan(b)]
print(b)
t = np.linspace(0, len(b), len(b)) #hours

np.savetxt('data_@1h_len{}.dat'.format(len(b)),b)
quit()

plt.plot(t, b)
plt.show()

quit()
a = []
for i in range(2):
	print(float(ds.data_vars['t2m'].time[i]))
	for k in range(data.shape[1]):
		t = float(ds.data_vars['t2m'][i,k,50,150])
		print(float(ds.data_vars['t2m'].step[k])/1e9, t)
		a.append(t)

plt.plot(a)
plt.show()
quit()
	#doing a latitude map
	#For some reasons there are a lot of nans around: probably this is required because not every resolution has the same angles??

lat_list = []
k = 2
for i in range(10,int(data.shape[2]/2)):
	long_list = []
	for j in range(10,int(data.shape[3]/2)):
		if not np.isnan(float(ds.data_vars['t2m'][1,k,i,j])):
			#print(k,i,j)
			#print(float(ds.data_vars['t2m'][1,k,i,j]))
			long_list.append(float(ds.data_vars['t2m'][1,k,i,j]))
	print(len(long_list))
	if len(long_list)==90:
		lat_list.append(long_list)


print(lat_list[0])
data_np = np.array(lat_list)
print(data_np.shape)
plt.imshow(data_np)
plt.show()


