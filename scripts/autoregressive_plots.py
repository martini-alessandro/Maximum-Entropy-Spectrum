import numpy as np
import matplotlib.pyplot as plt
from style_sheet import init_plotting

import sys
sys.path.insert(0,'..')
import memspectrum


N_process = 100 #10 simulated autoregressive processes
N_data = 50000 #number of timesteps for each autoregressive process
p_max = 5000

save_dir = 'arp_data/'
plot_dir = '../paper/Images/arp_errors/'

load = True
plot = True

	#initializing the models
true_arp_list = []
rec_MESA_list = []
data_list = []
loss_functions = ['FPE', 'CAT', 'OBD']

for i in range(N_process):

	if load:
		true_arp_list.append(memspectrum.MESA(save_dir+'arp_model_{}'.format(i)))
		#data_list.append(np.loadtxt(save_dir+'data_{}'.format(i)))
		rec_MESA_list.append( (memspectrum.MESA(save_dir+'MESA_model_{}_{}'.format(loss_functions[0],i)), 
			memspectrum.MESA(save_dir+'MESA_model_{}_{}'.format(loss_functions[1],i)),
			memspectrum.MESA(save_dir+'MESA_model_{}_{}'.format(loss_functions[2],i)) )
			)
	else:
		true_arp_list.append(memspectrum.MESA())
			#getting p
		p = np.random.choice(np.linspace(np.log10(2),np.log10(p_max), 10000))
		p = int(10**(p))
		
		a_k = np.random.normal(0, 0.1, size = (p,))
		a_k = np.multiply(a_k,np.exp(-0.05*np.arange(0,len(a_k),1)))
		a_k[0] = 1 #you forgot this part!!! Now you need to start again :(
		P = np.random.exponential(0.05)

		print('Process {} has p = {} and P, a_k = {}, {}'.format(i,p,P,None))
		
			#initializing the model
		true_arp_list[-1].N = N_data
		true_arp_list[-1].P = P
		true_arp_list[-1].a_k = a_k
		
		true_arp_list[-1].save(save_dir+'arp_model_{}'.format(i))
		
			#forecasting
		data_0 = np.random.rand(p) #initial data
		burn_steps = 1000*p
		
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
if plot:
	N_p_plot = 1 #N of series to be shown in the err vs p plot
	colors = {'FPE': 'r', 'CAT': 'k', 'OBD': 'green'}
	fig_p = init_plotting()
	ax_p = fig_p.gca()
	
	fig_scatter = init_plotting()
	ax_scatter = fig_scatter.gca()
	
		#computing erorrs
	p_rec = np.zeros((N_process,3))
	p_true = np.zeros((N_process,))
	diff_ak = np.zeros((N_process,N_data,3))
	for i in range(N_process):
		for j, l in enumerate(loss_functions):
			p_true[i] = true_arp_list[i].get_p()
			p_rec[i,j] = rec_MESA_list[i][j].get_p()#- true_arp_list[i].get_p() + rec_MESA_list[i][j].get_p()
			#print(j, true_arp_list[i].get_p() , rec_MESA_list[i][j].get_p(),p_rec[i,j])
			ak_true = np.concatenate([true_arp_list[i].a_k, np.zeros((N_data-len(true_arp_list[i].a_k),))]) #(N_data,)
			ak_rec = np.concatenate([rec_MESA_list[i][j].a_k, np.zeros((N_data-len(rec_MESA_list[i][j].a_k),))]) #(N_data,)
			diff_ak[i,:,j] = -(ak_true-ak_rec)
			

		#plotting
	for i, l in enumerate(loss_functions):
#		ax_scatter.scatter(p_rec[:,i], np.sqrt(np.mean(np.square(diff_ak[...,i]), axis = 1)), c = colors[l], label = l)
		ax_scatter.scatter(p_rec[:,i], p_true, c = colors[l], label = l, s = 2)
		ax_p.plot(range(N_data), diff_ak[:N_p_plot,:,i].T, 'o', c = colors[l])
		#ax_p.set_yscale('log')
	ax_scatter.legend(loc = 'upper left')
	ax_p.legend()
	ax_p.set_xlim([0,200])
	ax_p.set_ylim([-0.15,0.15])
	ax_p.set_xlabel("p")
	ax_p.set_ylabel("Difference in a_k")
	#ax_scatter.set_xlabel("p - p_true")
	#ax_scatter.set_ylabel("Averaged squared error")
	ax_scatter.set_xlabel(r"$p_{MESA}$")
	ax_scatter.set_ylabel(r"$p_{true}$")
	ax_scatter.plot(range(2,p_max),range(2,p_max))
	ax_scatter.set_xscale('log')
	ax_scatter.set_yscale('log')
	
	fig_p.tight_layout()
	fig_scatter.tight_layout()
	
	fig_scatter.savefig(plot_dir+'scatter_deltap_ptrue.pdf')
	fig_p.savefig(plot_dir+'scatter_deltaak_p.pdf')
	plt.show()




		
		
		
		
		
		
		
