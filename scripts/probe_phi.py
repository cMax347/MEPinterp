import numpy as np
import datetime
import os
from mep_worker import MEP_worker 
import matplotlib.pyplot as plt





class Phi_probe:

	def __init__(		self,n_phi, val_bands, mp_grid, gamma_scale, kubo_tol, 	
						hw, n_eF, eF_min, eF_max, Tkelvin, eta_smearing,debug_mode, do_gauge_trafo='T' ,
						do_write_velo='F', do_write_mep_bands='F',
						do_mep='T', do_kubo='F', do_ahc='F', do_opt='F', do_gyro='F'
						):
		self.n_phi 			= 	n_phi
		self.val_bands		= 	val_bands
		self.mp_grid		= 	mp_grid
		self.gamma_scale	=	gamma_scale
		self.kubo_tol		= 	kubo_tol
		self.hw				= 	hw 
		self.n_eF			=	n_eF
		self.eF_min			=	eF_min
		self.eF_max			=	eF_max
		self.Tkelvin		= 	Tkelvin 
		self.eta_smearing	= 	eta_smearing
		self.debug_mode		= 	debug_mode 
		self.do_gauge_trafo	= 	do_gauge_trafo
		self.do_write_velo	=	do_write_velo
		self.do_write_mep_bands	=	do_write_mep_bands
		self.do_mep			=	do_mep
		self.do_kubo		=	do_kubo
		self.do_ahc			=	do_ahc
		self.do_opt			=	do_opt
		self.do_gyro		=	do_gyro

		#derived attributes
		self.root_dir	= os.getcwd()+'/'+datetime.date.today().strftime("%d%B%Y")+'_mp'+str(self.mp_grid[0])+str(self.mp_grid[1])+str(self.mp_grid[2])
		
		#search new root folder name
		if os.path.isdir(self.root_dir):
			old_path 	= self.root_dir
			cnt 		= 0
			while os.path.isdir(old_path) and cnt < 5:
				old_path	= old_path + '.old'
				cnt			= cnt + 1
				try:
					os.rename(self.root_dir, old_path)
				except OSError:
					print(old_path+ ' exists already')

		self.plot_dir		= 	self.root_dir+'/mep_tot_plots'
		self.phi_tot_data 	= 	[]
		self.phi_cs_data 	= 	[]
		self.phi_lc_data 	= 	[]
		self.phi_ic_data 	= 	[]
		self.phi_bands_data	=	[]




	def __del__(self):
		print('done with Phi probing, by\n\n')






	def iterate_phi(self,mpi_np=1, plot_bandstruct=False):
		for phi  in np.linspace(0.0, 2.0, num = self.n_phi):		#iterate over relative phi (phi_rel = phi / np.pi)
			work_dir =	self.root_dir+'/phi'+str(phi)		
			phi_pi	 = 	phi
			#init current job
			worker = MEP_worker(	self.root_dir, 
									work_dir, 
									phi_pi, 						#give  phi_rel*pi to calculation
									self.val_bands, 
									self.mp_grid, 
									self.kubo_tol,
									self.hw,
									self.n_eF,
									self.eF_min,
									self.eF_max,
									self.Tkelvin,
									self.eta_smearing,
									self.debug_mode,
									self.do_gauge_trafo,
									self.do_write_velo,
									self.do_write_mep_bands,
									self.do_mep,
									self.do_kubo,
									self.do_ahc,
									self.do_opt,
									self.do_gyro
								)
			#run calc 
			worker.run(mpi_np=mpi_np)
			mep_tens, mep_cs, mep_lc, mep_ic, mep_bands = worker.get_mep_tens()
			self.phi_tot_data.append(		[phi, mep_tens	]			)
			self.phi_cs_data.append(		[phi, mep_cs	]			)
			self.phi_lc_data.append(		[phi, mep_lc	]			)
			self.phi_ic_data.append(		[phi, mep_ic	]			)
			self.phi_bands_data.append(		[phi, mep_bands ]			)
			#plot bands
			if plot_bandstruct:
				try:
					worker.plot_bands()
				except:
					print('could not plot bands')
			
		self.phi_tot_data 	= 	sorted(self.phi_tot_data)
		self.phi_cs_data	=	sorted(self.phi_cs_data)
		self.phi_lc_data	=	sorted(self.phi_lc_data)
		self.phi_ic_data	=	sorted(self.phi_ic_data)
		self.phi_bands_data	=	sorted(self.phi_bands_data)



	def print_results_container(self):
		print('^^^^^raw data ^^^^^^^')
		print('mep_bands :')
		print(self.phi_bands_data)
		#print(self.phi_bands_data[0][0])

		print('phi_tot_data =',self.phi_tot_data)
		print('------------------------------------')




	

	



				

	








def probe_phi(		n_phi, val_bands, mp_grid,mpi_np=1, gamma_scale=1,
					kubo_tol=1e-3, hw=0.001, n_eF=1,	eF_min=0.0, eF_max=0.0, Tkelvin=11.0, eta_smearing=0.2, 
					plot_bandstruct	=	False, 
					plot_orb_cont	=	True,
					plot_band_res	=	False,
					debug_mode		=	True, 
					do_gauge_trafo	=	True,
					do_write_velo	=	False,
					do_write_mep_bands	= False,
					do_mep			=	'T', 
					do_kubo			=	'F', 
					do_ahc			=	'F', 
					do_opt			=	'F', 
					do_gyro			=	'F'	
				):
	phi_job	= Phi_probe(	n_phi, val_bands, mp_grid, gamma_scale, 
							kubo_tol, hw, n_eF,	eF_min, eF_max, Tkelvin, eta_smearing,
							debug_mode, do_gauge_trafo, do_write_velo, do_write_mep_bands,
							do_mep, do_kubo, do_ahc, do_opt, do_gyro	
						)
	#
	phi_job.iterate_phi(plot_bandstruct=plot_bandstruct, mpi_np=mpi_np)
	phi_job.print_results_container()
	print("all done, by from probe_phi!")

	






probe_phi(		#parameter space probe density:
				n_phi				=	21					, 
				val_bands			=	2					, 
				#numerical parameters
				mp_grid				=	[16,16,16]			, 
				mpi_np				=	4					,
				gamma_scale			=	1.0					,
				kubo_tol			=	1e-5				, 
				hw					=	0.0					, 
				n_eF				=	1					,
				eF_min 				=	0.0					,
				eF_max				=	0.0					,
				Tkelvin				=	10.0				,
				eta_smearing		=	0.1					, 
				#set properties of plots
				plot_bandstruct		=	False				, 
				plot_orb_cont		=	True				,				
				plot_band_res		=	False				,
				#additional fortran controllers
				debug_mode			=	True				, 
				do_gauge_trafo		=	True				,	
				do_write_velo		=	False				,
				do_write_mep_bands	=	True				,
				do_mep				=	'T'					, 
				do_kubo				=	'F'					, 
				do_ahc				=	'F'					, 
				do_opt				=	'F'					, 
				do_gyro				=	'F'
			)



