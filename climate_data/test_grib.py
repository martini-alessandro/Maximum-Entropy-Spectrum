# https://github.com/ecmwf/cfgrib
# TO download data: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

ds = xr.open_dataset('adaptor.mars.internal-1614589133.4434555-10786-9-72b941a9-4c27-4d17-8061-3525cfcbd646.grib', engine='cfgrib', backend_kwargs={'read_keys': ['experimentVersionNumber']}) #20091231-20201130
#ds = xr.open_dataset('adaptor.mars.internal-1614619196.5842516-30992-8-350d6fd6-6cf6-4e7e-bed2-b92cca5f47ec.grib', engine='cfgrib', backend_kwargs={'read_keys': ['experimentVersionNumber']}) #19991231-20091231
#ds = xr.open_dataset('adaptor.mars.internal-1614619254.7736595-26399-6-0b5cc738-4574-4489-9cfe-8b807abaf151.grib', engine='cfgrib', backend_kwargs={'read_keys': ['experimentVersionNumber']}) #19891231-19991231


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

np.savetxt('data_@1h_len{}_.dat'.format(len(b)),b)
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


"""
array(['1999-12-31T00:00:00.000000000', '2000-01-01T00:00:00.000000000',
       '2000-01-02T00:00:00.000000000', ..., '2009-12-29T00:00:00.000000000',
       '2009-12-30T00:00:00.000000000', '2009-12-31T00:00:00.000000000'],
      dtype='datetime64[ns]')
Coordinates:
    number   int64 ...
  * time     (time) datetime64[ns] 1999-12-31 2000-01-01 ... 2009-12-31
    surface  int64 ...
Attributes:
    long_name:      initial time of forecast
    standard_name:  forecast_reference_time
<xarray.DataArray 'step' (step: 24)>
array([ 3600000000000,  7200000000000, 10800000000000, 14400000000000,
       18000000000000, 21600000000000, 25200000000000, 28800000000000,
       32400000000000, 36000000000000, 39600000000000, 43200000000000,
       46800000000000, 50400000000000, 54000000000000, 57600000000000,
       61200000000000, 64800000000000, 68400000000000, 72000000000000,
       75600000000000, 79200000000000, 82800000000000, 86400000000000],
      dtype='timedelta64[ns]')
Coordinates:
    number   int64 ...
  * step     (step) timedelta64[ns] 01:00:00 02:00:00 ... 1 days 00:00:00
    surface  int64 ...
Attributes:
    long_name:      time since forecast_reference_time
    standard_name:  forecast_period
<xarray.DataArray 'latitude' (latitude: 2)>
array([45.5, 45.4])
Coordinates:
    number    int64 ...
    surface   int64 ...
  * latitude  (latitude) float64 45.5 45.4
Attributes:
    units:             degrees_north
    standard_name:     latitude
    long_name:         latitude
    stored_direction:  decreasing
<xarray.DataArray 'longitude' (longitude: 2)>
array([9.1, 9.2])
Coordinates:
    number     int64 ...
    surface    int64 ...
  * longitude  (longitude) float64 9.1 9.2
Attributes:
    units:          degrees_east
    standard_name:  longitude
    long_name:      longitude
(3654, 24, 2, 2)
(3654, 24)
[[      nan       nan       nan ...       nan       nan 270.97482]
 [270.26868 269.6307  269.4444  ... 270.64658 270.48184 270.25766]
 [269.94144 269.96716 270.04718 ... 269.5512  269.8716  269.7848 ]
 ...
 [272.35785 271.99747 271.58832 ... 273.0402  273.2572  273.31464]
 [273.13177 272.9658  272.70944 ... 272.12778 271.6813  271.66272]
 [271.39862 271.34235 271.3499  ... 275.75543 275.4552        nan]]
[      nan       nan       nan ... 275.75543 275.4552        nan]
[270.97482 270.26868 269.6307  ... 275.48245 275.75543 275.4552 ]
"""


