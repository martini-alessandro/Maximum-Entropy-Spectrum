import sys
sys.path.insert(0,'..')
import mesa
import matplotlib.pyplot as plt
import numpy as np

import mlgw.GW_generator as gen
t = np.linspace(-10.,0.1,4000)
g = gen.GW_generator()
h_p, h_c = g.get_WF([20,10,0.1,-0.5],t)

#data = h_p[:200]*1e20
data = np.sin(10.*t)

#data = np.random.normal(0,1,10000)*np.random.poisson(10,10000) #this is a sanity check that it works

plt.plot(data)
plt.show()
m = mesa.MESA(data)

m.solve(method = "Fast", optimisation_method = "FPE")

np.random.seed(0)
old_f = m.forecast(500,1)

np.random.seed(0)
vec_f = m.forecast_vectorized(500,1, include_data=False)

print(np.allclose(vec_f,old_f))

plt.plot(np.squeeze(vec_f), label = 'vec')
plt.plot(np.squeeze(old_f), label = 'old')
plt.legend()
plt.show()
