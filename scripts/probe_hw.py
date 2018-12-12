import numpy as np
import datetime
import os
from mep_worker import MEP_worker 
import matplotlib.pyplot as plt





class hw_probe:




	def __init__(		self,
						tb_model,
						hw_min,
						hw_max,
						n_hw, 
						val_bands, 
						mp_grid, 
						gamma_scale, 
						kubo_tol, 	
						laser_phase, 
						eFermi, Tkelvin, eta_smearing,
						debug_mode='F', 
						do_gauge_trafo='T',
						R_vect_float='F',
						do_write_velo='F',
						do_mep='T', do_kubo='F', do_ahc='F', do_opt='F', do_gyro='F'
						):
		self.tb_model		=	tb_model
		self.hw_min			=	hw_min
		self.hw_max			=	hw_max
		self.n_hw 			= 	n_hw
		self.val_bands		= 	val_bands
		self.mp_grid		= 	mp_grid
		self.gamma_scale	=	gamma_scale
		self.kubo_tol		= 	kubo_tol
		self.laser_phase	=	laser_phase
		self.eFermi			= 	eFermi 
		self.Tkelvin		= 	Tkelvin 
		self.eta_smearing	= 	eta_smearing
		self.debug_mode		= 	debug_mode 
		self.do_gauge_trafo	= 	do_gauge_trafo
		self.R_vect_float	=	R_vect_float
		self.do_write_velo	=	do_write_velo
		self.do_mep			=	do_mep
		self.do_kubo		=	do_kubo
		self.do_ahc			=	do_ahc
		self.do_opt			=	do_opt
		self.do_gyro		=	do_gyro
		self.t_start		=	datetime.date.today()

		print("start hw-probe at "+self.t_start.strftime("%d%B%Y"))
		#derived attributes
		self.root_dir	= os.getcwd()+'/'+self.t_start.strftime("%d%B%Y")+'_mp'+str(self.mp_grid[0])+str(self.mp_grid[1])+str(self.mp_grid[2])
		
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
		#
		self.phi_tot_data = []
		self.phi_cs_data = []
		self.phi_lc_data = []
		self.phi_ic_data = []




	def __del__(self):
		time_spent	=	datetime.date.today() - self.t_start
		print('done with hw probing after '+str(time_spent)+', by\n\n')





	def iterate_hw(self,mpi_np=1, plot_bands=False):
		for hw  in np.linspace(self.hw_min, self.hw_max, num = self.n_hw):		#iterate over relative phi (phi_rel = phi / np.pi)
			work_dir =	self.root_dir+'/hw'+str(hw)		
			phi_pi	 = 	0.0

			print(mpi_np)
			#init current job
			# mep_worker  constructor
	#def __init__(		self, tb_model, root_dir, work_dir, phi, val_bands, mp_grid, 
	#					kubo_tol=1e-3,  hw=0.0, eFermi=0.0, Tkelvin=0.0, eta_smearing=0.0, 
	#					debug_mode='F', do_gauge_trafo='T',	
	#					do_write_velo='F', do_write_mep_bands='F',
	#					do_mep='T', do_kubo='F', do_ahc='F', do_opt='F', do_gyro='F'
	#					):



			worker = MEP_worker(	self.tb_model,
									self.root_dir, 
									work_dir, 
									phi_pi, 						#give  phi_rel*pi to calculation
									self.val_bands, 
									self.mp_grid, 
									self.kubo_tol,
									hw,
									self.laser_phase,
									self.eFermi,
									self.Tkelvin,
									self.eta_smearing,
									self.debug_mode,
									self.do_gauge_trafo,
									self.R_vect_float,
									self.do_write_velo,
									self.do_mep,
									self.do_kubo,
									self.do_ahc,
									self.do_opt,
									self.do_gyro
								)
			#run calc 
			worker.run(mpi_np=mpi_np)
			mep_tens, mep_cs, mep_lc, mep_ic, mep_bands 	= 	worker.get_mep_tens()
			ahc_tens										=	worker.get_ahc_tens()
			ohc_tens										=	worker.get_ohc_tens()
			
			self.phi_tot_data.append(		[hw, mep_tens	]			)
			self.phi_cs_data.append(		[hw, mep_cs	]			)
			self.phi_lc_data.append(		[hw, mep_lc	]			)
			self.phi_ic_data.append(		[hw, mep_ic	]			)
			#plot bands
			if plot_bands:
				try:
					worker.plot_bands()
				except:
					print('could not plot bands')
			
		self.phi_tot_data 	= 	sorted(self.phi_tot_data)
		self.phi_cs_data	=	sorted(self.phi_cs_data)
		self.phi_lc_data	=	sorted(self.phi_lc_data)
		self.phi_ic_data	=	sorted(self.phi_ic_data)



	def print_results_container(self):
		print('results =',self.phi_tot_data)










def probe_hw(		tb_model		=	"FeMn3q" 		,
					hw_min			=	0				,
					hw_max			= 	6				,
					n_hw			= 	10				,
					val_bands		=	2				, 
					mp_grid			= 	[1,1,1]			, 
					mpi_np			=	1				, 
					gamma_scale		=	1				,
					kubo_tol		=	1e-3			, 
					laser_phase		=	2.0				,
					eFermi			=	0.0				, 
					Tkelvin			=	11.0			, 
					eta_smearing	=	0.2				, 
					plot_bands		=	False			, 
					debug_mode		=	True			, 
					do_gauge_trafo	=	True			,
					R_vect_float	=	False			,
					do_write_velo	=	False			,
					do_mep			=	'T'				, 
					do_kubo			=	'F'				, 
					do_ahc			=	'F'				, 
					do_opt			=	'F'				, 
					do_gyro			=	'F'					
				):
	myTest	= hw_probe(	tb_model,hw_min, hw_max, n_hw, val_bands, mp_grid, gamma_scale, 
							kubo_tol, laser_phase, eFermi, Tkelvin, eta_smearing,
							debug_mode, do_gauge_trafo, R_vect_float,  do_write_velo,
							do_mep, do_kubo, do_ahc, do_opt, do_gyro	
						)
	#
	# mep_worker  constructor
	#def __init__(		self, tb_model, root_dir, work_dir, phi, val_bands, mp_grid, 
	#					kubo_tol=1e-3,  hw=0.0, eFermi=0.0, Tkelvin=0.0, eta_smearing=0.0, 
	#					debug_mode='F', do_gauge_trafo='T',	
	#					do_write_velo='F', do_write_mep_bands='F',
	#					do_mep='T', do_kubo='F', do_ahc='F', do_opt='F', do_gyro='F'
	#					):


	myTest.iterate_hw( mpi_np=mpi_np, plot_bands=plot_bands)
	myTest.print_results_container()



probe_hw(		tb_model		= "FeMn3q"			,
				hw_min			= 0					,
				hw_max			= 6					,
				n_hw			= 7					,
				val_bands		= 2					, 
				mp_grid			= [16,16,16]		, 
				mpi_np			= 4					,
				gamma_scale		= 7.7481e-5			,
				kubo_tol		= 1e-5				, 
				laser_phase		= 1.0				,
				eFermi			= 0.0				, 
				Tkelvin			= 300.0				, 
				eta_smearing	= 0.1				, 
				debug_mode		= True				, 
				do_gauge_trafo	= True				,
				R_vect_float	= True				,	
				do_write_velo	= False				,
				do_mep			= 'T'				, 
				do_kubo			= 'F'				, 
				do_ahc			= 'F'				, 
				do_opt			= 'F'				, 
				do_gyro			= 'F'
			)