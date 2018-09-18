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


	def add_jobs(self, phi, val_bands, mp_dens_per_dim, hw, eFermi, Tkelvin, eta_smearing):
		for n_mp in mp_dens_per_dim:
			nK 		=	n_mp**3
			mp_grid	=	[n_mp, n_mp, n_mp]

			work_dir= self.root_dir+'/nK'+str(nK)
			self.work_dirs.append(work_dir)

			job = MEP_worker(self.root_dir, work_dir, phi, val_bands, mp_grid, hw, eFermi, Tkelvin, eta_smearing	)
			self.jobs.append( 	job	)


	def run_jobs(self, mpi_np=1):
		for job in self.jobs:
			job.run(mpi_np)



#************************************************************************************************************************************************************************





def test(hw, eFermi, Tkelvin, eta_smearing):
	root_dir	=	os.getcwd()+'/'+datetime.date.today().strftime("%d%B%Y")+'_k_conv_run_TEST'

	test	= conv_run(root_dir)
	mp_dens = [1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 80 ,96, 128]
	test.add_jobs(0.0, 2, mp_dens, hw, eFermi, Tkelvin, eta_smearing)

	test.run_jobs(mpi_np=16)




def run_cluster(hw, eFermi, Tkelvin, eta_smearing):
	root_dir		=	os.getcwd()+'/k_conv_cluster'


	val_bands		=	2
	mp_dens			=	[1, 2, 4, 6, 8, 12, 16]#, 42, 32, 48, 64, 80, 128,256]
	phi_lst			=	[0.0, .1]

	for phi in phi_lst:
		cluster_calc 	= 	conv_run(root_dir+'_phi'+str(phi))
		#
		cluster_calc.add_jobs(phi,	val_bands,	mp_dens, hw, eFermi, Tkelvin, eta_smearing)
		cluster_calc.run_jobs(mpi_np=16)





#************************************************************************************************************************************************************************

#do_test	=	input("do you want to run a test job or prepare cluster? (t/c)")
#breaker	= 0
#while (breaker<10) 	and 	(do_test is not "t")  and (do_test is not "c"): 
#	print('input has to be "t" or "c"')
#	breaker = breaker+ 1
#	do_test	=	input("do you want to run a test job or prepare cluster? (t/c)")


#if do_test is "t":
#	test()
#elif do_test is "c":
#	prepare_cluster()
#else:
#	print("I dont know what I should do, so I do nothing instead")

hw				= 	0.0
eFermi			=	0.0	
Tkelvin			=	0.0	
eta_smearing	=	0.0



run_cluster(hw,eFermi, Tkelvin, eta_smearing)


