import numpy as np
import datetime
import os
from mep_worker import MEP_worker 
import matplotlib.pyplot as plt



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#							CLASS TO HANDLE HW (LASERFIELD) PARAMETERSPACE PROBES
class hw_job:
#
#		------------------------------------------------------------------------------
#
#
	def __init__(		self, tb_model, use_pos_op,  phi_pi, val_bands, mp_grid, gamma_scale, kubo_tol, 	
						laser_phase, n_eF, eF_min, eF_max, Tkelvin, eta_smearing,debug_mode, do_gauge_trafo='T' ,
						do_r_vect_float='T',
						do_write_velo='F', do_write_mep_bands='F',
						do_mep='T', do_kubo='F', do_ahc='F', do_opt='F', do_gyro='F'
						):
		self.tb_model		=	tb_model
		self.use_pos_op		=	use_pos_op
		self.phi_pi			=	phi_pi
		self.val_bands		= 	val_bands
		self.mp_grid		= 	mp_grid
		self.gamma_scale	=	gamma_scale
		self.kubo_tol		= 	kubo_tol
		self.laser_phase	=	laser_phase
		self.n_eF			=	n_eF
		self.eF_min			=	eF_min
		self.eF_max			=	eF_max
		self.Tkelvin		= 	Tkelvin 
		self.eta_smearing	= 	eta_smearing
		self.debug_mode		= 	debug_mode 
		self.do_gauge_trafo	= 	do_gauge_trafo
		self.do_r_vect_float=	do_r_vect_float
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
#
#	---------------
#
	def __del__(self):
		print('done with hw probing, by\n\n')
#
#	---------------
#

	def probe_hw_space(self,hw_min, hw_max, n_hw, mpi_np=1, plot_bandstruct=False):
		print('[probe_hw_space]: hw_min	= '+str(hw_min)+'	hw_max='+str(hw_max)+'	n_hw='+str(n_hw)	)
		for hw  in np.linspace(hw_min, hw_max, num = n_hw):		#iterate over relative phi (phi_rel = phi / np.pi)
			work_dir =	self.root_dir+'/hw'+str(hw)		
			#init current job
			print("")
			print("")
			print("")
			print("")
			print('[probe_hw_space]: init new MEP_worker (hw='+str(hw)+')')
			worker = MEP_worker(	self.tb_model,
						self.use_pos_op,
						self.root_dir,
						work_dir, 
						self.phi_pi, 						#give  phi_rel*pi to calculation
						self.val_bands, 
						self.mp_grid, 
						self.kubo_tol,
						hw,
						self.laser_phase,
						self.n_eF,
						self.eF_min,
						self.eF_max,
						self.Tkelvin,
						self.eta_smearing,
						self.debug_mode,
						self.do_gauge_trafo,
						self.do_r_vect_float,
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
			print('[probe_hw_space]: curr. calculation finished, try to collect data ...')
			mep_tens, mep_cs, mep_lc, mep_ic, mep_bands = worker.get_mep_tens()
			self.phi_tot_data.append(		[hw, mep_tens	]			)
			self.phi_cs_data.append(		[hw, mep_cs	]			)
			self.phi_lc_data.append(		[hw, mep_lc	]			)
			self.phi_ic_data.append(		[hw, mep_ic	]			)
			self.phi_bands_data.append(		[hw, mep_bands ]			)
			#plot bands
			if plot_bandstruct:
				try:
					worker.plot_bands()
				except:
					print('[probe_hw_space]: could not plot bands')
			print('[probe_hw_space]: finished  hw='+str(hw)+' MEP_worker')
		#
		self.phi_tot_data 	= 	sorted(self.phi_tot_data)
		self.phi_cs_data	=	sorted(self.phi_cs_data)
		self.phi_lc_data	=	sorted(self.phi_lc_data)
		self.phi_ic_data	=	sorted(self.phi_ic_data)
		self.phi_bands_data	=	sorted(self.phi_bands_data)
		print("*")
		print("*")
		print("*")

#
#	---------------
#
	def print_results_container(self):
		print('^^^^^raw data ^^^^^^^')
		print('mep_bands :')
		print(self.phi_bands_data)
		#print(self.phi_bands_data[0][0])

		print('phi_tot_data =',self.phi_tot_data)
		print('------------------------------------')
#
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#---------------------------------------------------------------------------------------------------------------------------------------------------------
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




	

	



				

	
#	
#~~~~~~~~~		CONTROL TB-PARAMETERS BELLOW	~~~~~~~~~~~~~~~~~
#
#	I. setup new model system
tb_system = hw_job(		#parameter space probe density:
								tb_model			=	'FeMn3q'			,
								use_pos_op			=	False				,
								#sys para
								phi_pi				=	0					,
								val_bands			=	2					, 
								#numerical parameters
								mp_grid				=	[16,16,16]			, 
								gamma_scale			=	1.0					,
								kubo_tol			=	1e-5				, 
								laser_phase			=	1.0					,
								n_eF				=	1					,
								eF_min 				=	0.0					,
								eF_max				=	0.0					,
								Tkelvin				=	10.0				,
								eta_smearing		=	0.1					, 
								#additional fortran controllers
								debug_mode			=	False				, 
								do_gauge_trafo		=	True				,	
								do_r_vect_float		=	True				,
								do_write_velo		=	False				,
								do_write_mep_bands	=	True				,
								do_mep				=	'F'					, 
								do_kubo				=	'F'					, 
								do_ahc				=	'T'					, 
								do_opt				=	'T'					, 
								do_gyro				=	'T'
					)				
#
#	II.	probe hw between hw_min & hw_max
tb_system.probe_hw_space(		hw_min			=		0		, 
								hw_max			=		6		,	
								n_hw			= 		7		,
								plot_bandstruct =   	True	, 
								mpi_np			=		4
					)
#
#	III. show results in std out
tb_system.print_results_container()
#
print("all done, by from probe_hw!")










