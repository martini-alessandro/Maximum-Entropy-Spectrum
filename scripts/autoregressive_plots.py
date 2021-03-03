import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'..')
import memspectrum



N_process = 10 #10 simulated autoregressive processes
N_data = 1000 #number of timesteps for each autoregressive process
p_max = 100

save_dir = 'arp_data/'

srate = 1.

load = True

	#initializing the models
true_arp_list = []
rec_MESA_list = []
data_list = []
loss_functions = ['FPE', 'CAT', 'OBD']

for i in range(N_process):

	if load:
		true_arp_list.append(memspectrum.MESA(save_dir+'arp_model_{}'.format(i)))
		data_list.append(np.loadtxt(save_dir+'data_{}'.format(i)))
		rec_MESA_list.append( (memspectrum.MESA(save_dir+'MESA_model_{}_{}'.format(loss_functions[0],i)), 
			memspectrum.MESA(save_dir+'MESA_model_{}_{}'.format(loss_functions[1],i)),
			memspectrum.MESA(save_dir+'MESA_model_{}_{}'.format(loss_functions[2],i)) )
			)
	else:
		true_arp_list.append(memspectrum.MESA())
			#getting p
		p = np.random.choice(range(2,p_max))
		
		a_k = np.random.normal(0, 0.2, size = (p,))
		a_k = np.multiply(a_k,np.exp(-0.1*np.arange(0,len(a_k),1)))
		P = np.random.exponential(0.05)

		print('Process {} has p = {} and P, a_k = {}, {}'.format(i,p,P,a_k))
		
			#initializing the model
		true_arp_list[-1].N = N_data
		true_arp_list[-1].P = P
		true_arp_list[-1].a_k = a_k
		
		true_arp_list[-1].save(save_dir+'arp_model_{}'.format(i))
		
			#forecasting
		data_0 = np.random.rand(p) #initial data
		burn_steps = 100000
		
		data = true_arp_list[-1].forecast(data_0, burn_steps+N_data)[0,:] #(burn_steps+N_data,)
		
		data_list.append(data[-N_data:])
		#print("Time series: ", data_list[-1])
		
		#plt.plot(data_list[-1])
		#plt.show()
		np.savetxt(save_dir+'data_{}'.format(i), data_list[-1])
		
		
			#reconstructing with MESA
		temp_MESA =  (memspectrum.MESA(), memspectrum.MESA(), memspectrum.MESA()) 
		for j, l in enumerate(loss_functions):
			print("Solving with {} @ iteration {}".format(l,i))
			temp_MESA[j].solve(data_list[-1], method = 'standard', early_stop = False, optimisation_method = l)
			temp_MESA[j].save(save_dir+'MESA_model_{}_{}'.format(l,i))
		rec_MESA_list.append(temp_MESA)



	#doing plots






		
		
		
		
		
		
		
