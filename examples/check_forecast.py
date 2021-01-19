try:
	import sys
	sys.path.insert(0,'..')
	import mesa
except:
	import mesa

import matplotlib.pyplot as plt
import numpy as np

import mlgw.GW_generator as gen
t = np.linspace(-10.,0.1,4000)
g = gen.GW_generator()
h_p, h_c = g.get_WF([20,10,0.1,-0.5],t)

data = h_p[-3000:-2000]*1e20
#data = np.sin(10.*t)

#data = np.random.normal(0,1,10000)*np.random.poisson(10,10000) #this is a sanity check that it works

plt.figure()
plt.plot(data)
#plt.show()

m = mesa.MESA(data)

m.solve(method = "Fast", optimisation_method = "FPE")

np.random.seed(0)
forecast = m.forecast(1500,1, include_data = True)

print(len(m.a_k))
plt.figure()
plt.plot(np.squeeze(forecast[0,:]), label = 'vec')
plt.legend()


plt.figure()
plt.plot(*m.spectrum(np.diff(t)[0]))


plt.show()
