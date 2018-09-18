import numpy as np
import datetime
import os
from mep_worker import MEP_worker 
import matplotlib.pyplot as plt


from fortran_io			import 	read_mep_file





class Phi_probe:

	def __init__(self,n_phi, val_bands, mp_grid ):
		self.n_phi 			= n_phi
		self.val_bands		= val_bands
		self.mp_grid		= mp_grid
		#derived attributes
		self.root_dir	= os.getcwd()+'/'+datetime.date.today().strftime("%d%B%Y")+'_mp'+str(self.mp_grid[0])+str(self.mp_grid[1])+str(self.mp_grid[2])
		
	
		self.plot_dir	= self.root_dir+'/mep_tot_plots'
	




	def __del__(self):
		print('plotted Phi probing, by\n\n')



	def print_results_container(self):
		print('results =',self.phi_tot_data)




	def read_phi(self):
		self.phi_tot_data = []
		self.phi_cs_data = []
		self.phi_lc_data = []
		self.phi_ic_data = []
		#
		for dirName, subdirList, fileList in os.walk(self.root_dir):
			print("SUBDIRLIST="+str(subdirList))
			if 'phi' in dirName and not 'bands' in dirName:
				print('Found directory: %s' % dirName)
				print("associated phi= "+dirName.split("phi")[1])

		for phi  in np.linspace(0.0, 2.0, num = self.n_phi):		#iterate over relative phi (phi_rel = phi / np.pi)
			work_dir =	self.root_dir+'/phi'+str(phi)		
			phi_pi	 = 	phi
			
			mep_file_path	= work_dir+'/mep/mep_tens.dat'
			mep_tens		= read_mep_file(mep_file_path)
			#CHERN-SIMONS
			mep_file_path	= work_dir+'/mep/mep_cs.dat'
			mep_cs			= read_mep_file(mep_file_path)
			#LOCAL
			mep_file_path	= work_dir+'/mep/mep_lc.dat'
			mep_lc		= read_mep_file(mep_file_path)
			#ITINERANT
			mep_file_path	= work_dir+'/mep/mep_ic.dat'
			mep_ic		= read_mep_file(mep_file_path)

			#only record if the container is not empty (i.e. the file was found and had good behaviour)
			if len(mep_tens) is 3:
				self.phi_tot_data.append(		[phi, mep_tens	]			)
			else:
				print("skip mep_tens for phi="+str(phi))
				print("len(mep_tens)="+str(len(mep_tens)))
			if len(mep_cs) is 3:
				self.phi_cs_data.append(		[phi, mep_cs	]			)
			if len(mep_lc) is 3:
				self.phi_lc_data.append(		[phi, mep_lc	]			)
			if len(mep_ic) is 3:
				self.phi_ic_data.append(		[phi, mep_ic	]			)
		#			
		self.phi_tot_data 	= 	sorted(self.phi_tot_data)
		self.phi_cs_data	=	sorted(self.phi_cs_data)
		self.phi_lc_data	=	sorted(self.phi_lc_data)
		self.phi_ic_data	=	sorted(self.phi_ic_data)




	









	

	def plot_mep_over_phi(self, plot_contributions=True, label_size=14, xtick_size=12, ytick_size=12):
		try:
			os.mkdir(self.plot_dir)
		except OSError:
			print('Could not make directory ',self.plot_dir)

		phi_plot 		= []
		mep_tot_data 	= []
		mep_cs_data		= []
		mep_lc_data		= []
		mep_ic_data		= []
		for data in self.phi_tot_data:
			phi_plot.append(data[0])
			mep_tot_data.append(data[1])
		
		if plot_contributions:
			for data in self.phi_cs_data:
				mep_cs_data.append(data[1])
			for data in self.phi_lc_data:
				mep_lc_data.append(data[1])
			for data in self.phi_ic_data:
				mep_ic_data.append(data[1])


		mep_max		= np.amax(mep_tot_data)
		mep_min 	= np.amin(mep_tot_data)
		mep_delta 	= (mep_max - mep_min) / 100.0
		mep_max		= mep_max + mep_delta
		mep_min		= mep_min - mep_delta




		dim_str	= []
		dim_str.append('x')
		dim_str.append('y')
		dim_str.append('z')

		#plot all 9 components of mep tensor
		for i in range(0,3):
			for j in range(0,3):
				#COLLECT a_ij(phi=0:n_phi)
				mep_tot_plot 	= []
				mep_cs_plot		= []
				mep_lc_plot		= []
				mep_ic_plot		= []
				for mep_tens in mep_tot_data:
					mep_tot_plot.append(	mep_tens[i][j]	)
				if plot_contributions:
					for mep_cs in mep_cs_data:
						mep_cs_plot.append(		mep_cs[i][j]	)
					for mep_lc in mep_lc_data:
						mep_lc_plot.append(		mep_lc[i][j]	)
					for mep_ic in mep_ic_data:
						mep_ic_plot.append(		mep_ic[i][j]	)






				#do PLOT
				fig, ax  = plt.subplots(1,1) 
				plt.plot(phi_plot, mep_tot_plot,'-+', color='black',label="tot")

				if plot_contributions:
					plt.plot(phi_plot, mep_cs_plot,		'--', 	color='red',		label="CS")
					plt.plot(phi_plot, mep_lc_plot,		'o-', 	color='darkblue',	label="LC")
					plt.plot(phi_plot, mep_ic_plot,		'^-', 	color='lightblue',	label="IC")


				
				#X-AXIS
				plt.xlabel(r'$\varphi$',	fontsize=label_size)
				ax.set_xlim([0,2])
				ax.set_xticks(np.array([0,0.5,1,1.5,2]))
				ax.set_xticklabels(np.array([r'$0$','',r'$\pi$','',r'$2\pi$']))
				plt.tick_params(axis='x',which='major', direction='in',labelsize=xtick_size)
				
				#Y-AXIS
				plt.ylabel(r'$\alpha_{'+dim_str[i]+dim_str[j]+'}$',	fontsize=label_size)
				ax.set_ylim([mep_min,mep_max])
				plt.tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)

				if plot_contributions:
					plt.legend()

				plt.tight_layout()
				plt.savefig(self.plot_dir+'/mep_'+dim_str[i]+dim_str[j]+'.pdf')
				plt.close()


				

	








def unit_test(n_phi, val_bands, mp_grid):
	myTest	= Phi_probe(n_phi, val_bands, mp_grid)
	#
	myTest.read_phi()
	myTest.print_results_container()
	
	myTest.plot_mep_over_phi(label_size=14, xtick_size=12, ytick_size=12)







unit_test(n_phi=11, val_bands=2, mp_grid=[128,128,128]	)



