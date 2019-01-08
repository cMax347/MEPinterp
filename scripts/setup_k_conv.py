import numpy as np
import datetime
import os

from mep_worker import MEP_worker




#************************************************************************************************************************************************************************

class conv_run:

	def __init__(self, root_dir):
		self.work_dirs	= []
		self.jobs		= []
		self.root_dir	= root_dir
		#
		#create root directory
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
		os.mkdir(self.root_dir)


	def add_jobs(	self, 	
					tb_model, use_pos_op, phi, val_bands, mp_dens_per_dim, 
					kubo_tol, hw, laser_phase,  n_eF, eF_min, eF_max,  Tkelvin, eta_smearing,
					debug_mode, do_gauge_trafo='T',	R_vect_float='F', do_write_velo='F',
					do_mep='T', do_kubo='F', do_ahc='F', do_opt='F', do_gyro='F'
				):
		for n_mp in mp_dens_per_dim:
			nK 		=	n_mp**3
			mp_grid	=	[n_mp, n_mp, n_mp]

			work_dir= self.root_dir+'/nK'+str(nK)
			self.work_dirs.append(work_dir)



			job = MEP_worker(	tb_model, use_pos_op, self.root_dir, work_dir, phi, val_bands, mp_grid, 
								kubo_tol, hw, laser_phase, n_eF, eF_min, eF_max, Tkelvin, eta_smearing, 
								debug_mode, do_gauge_trafo,	R_vect_float,	do_write_velo,
								do_mep, do_kubo, do_ahc, do_opt, do_gyro	
							)
			self.jobs.append( 	job	)


	def run_jobs(self, mpi_np=1):
		for job in self.jobs:
			job.run(mpi_np)



#************************************************************************************************************************************************************************
#************************************************************************************************************************************************************************
#************************************************************************************************************************************************************************





#paras
root_dir		=	os.getcwd()+'/k_conv_'+datetime.date.today().strftime("%d%B%Y")

#Laser
hw				= 	1.0
eta_smearing	=	0.1
laser_phase		=	1.0




#FERMI SMEARING
n_eF			=	25 
eF_min			=	-5.0
eF_max			=	-4.5	
Tkelvin			=	300.0	


#flags
debug_mode		=	'F'
do_gauge_trafo	=	'T'
do_write_velo	=	'F'

#repsonse tensors to calculate:
kubo_tol		=	1e-5
do_mep			=	'T'
do_kubo			=	'T'
do_ahc			=	'T'
do_opt			=	'T'
do_gyro			=	'F'

#TB MODEL
tb_model		=	'FeMn3q'
use_pos_op		=	False
R_vect_float	=	True


#NUMERICS
val_bands		=	1
mp_dens			=	[1, 2, 4]#, 6, 8, 12, 16]#, 32, 48, 64, 80, 128,256,512]
phi_lst			=	[0.0] 	#,1.0,2.0]
n_mpi_procs		=	4


conv_run	 	= 	conv_run(root_dir)
conv_run.add_jobs(	tb_model, use_pos_op, phi,	val_bands,	mp_dens, 
							kubo_tol, hw, laser_phase,  n_eF, eF_min, eF_max,  Tkelvin, eta_smearing, 
							debug_mode, do_gauge_trafo, R_vect_float, do_write_velo,
							do_mep, do_kubo, do_ahc, do_opt, do_gyro	
						)
conv_run.run_jobs(mpi_np=n_mpi_procs)

#for phi in phi_lst:
#	cluster_calc 	= 	conv_run(root_dir+'_phi'+str(phi))
#	#
#	cluster_calc.add_jobs(	tb_model, use_pos_op, phi,	val_bands,	mp_dens, 
#							kubo_tol, hw, laser_phase,  n_eF, eF_min, eF_max,  Tkelvin, eta_smearing, 
#							debug_mode, do_gauge_trafo, R_vect_float, do_write_velo,
#							do_mep, do_kubo, do_ahc, do_opt, do_gyro	
#						)
#	cluster_calc.run_jobs(mpi_np=n_mpi_procs)



