import numpy as np
import matplotlib.pyplot as plt
try:
	import sys
	sys.path.insert(0,'..')
	import mesa
	import GenerateTimeSeries
except:
	import mesa
	import mesa.GenerateTimeSeries as GenerateTimeSeries

N, dt = 1000, .01 		#the number of samples and sampling rate 
f = 2 	 	   	    	#The frequency of the sinusoidal signal 
time = np.arange(0, N) * dt 	#The time array 
data = np.sin(2 * np.pi * f * time) + np.random.normal(scale = 0.4,
                                                       size = 1000)

M = mesa.MESA(data) 		#Initialize MESA class
P, a_k, opt = M.solve() 	#Solve method and returns the 

f = np.linspace(0,50,10000)

M.save_spectrum("spectrum.dat", dt)

spec = np.loadtxt('spectrum.dat')
print(spec.shape)
plt.plot(spec[:,0], spec[:,1])
plt.show()
