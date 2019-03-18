import os
import time
import shutil as sh
import numpy as np
#
from scipy.optimize import 	minimize
from scipy.optimize import	basinhopping 	
#
from FeMn_WX_rashba import 	write_3q_HR
from fortran_io		import	read_real_tens_file


def save_mkdir(new_dir):
	if os.path.isdir('./'+new_dir):
		sh.rmtree(new_dir)	
	os.mkdir(new_dir)


def run_fortran(run_dir,	N_mpi_jobs=1,N_omp_jobs=1):
	start 		= time.time()
	root_dir	=	os.getcwd()
	os.chdir(run_dir)
	#
	os.environ['OMP_NUM_THREADS'] = str(N_omp_jobs)
	os.system('mpirun -np '+str(N_mpi_jobs)+' ./mepInterp'+'	> mep.log 2> mep.err ')
	#
	os.chdir(root_dir)
	end			=	time.time()
	print("[run_fortran]:	finished this run after ",end-start," sec")



def get_constraints():
	#	x:	t, strain, J_ex, tso, theta_deg
	constr = []
	#	hopping
	constr.append({'type':'ineq', 'fun' : lambda x: 	0.0-x[0]			})		#		-2 <= x[0]	<= 0.1
	constr.append({'type':'ineq', 'fun' : lambda x: 	2-np.abs(x[0])		})				 
	#	strain
	constr.append({'type':'ineq', 'fun' : lambda x: 	x[1]-.7 	})		#		.5	<= x[1] <= 1.5
	constr.append({'type':'ineq', 'fun' : lambda x: 	1.3- np.abs(x[1])	})
	#	J_ex
	constr.append({'type':'ineq', 'fun' : lambda x: 	-x[2]-2.0 	})		#		
	constr.append({'type':'ineq', 'fun' : lambda x: 	6.0 - np.abs(x[2])	})
	#	tso
	constr.append({'type':'ineq', 'fun' : lambda x: 	 - x[3]	})
	constr.append({'type':'ineq', 'fun' : lambda x: 	.6 - np.abs(x[3])	})
	#
	constr.append({'type':'ineq', 'fun' : lambda x: 	x[4] 		})
	constr.append({'type':'ineq', 'fun' : lambda x: 	90.1-np.abs(x[4])	})
	#
	return constr


class MEP_optimizer:

	def __init__(	self, 
					N_mpi_procs=1,
					base_dir='./mep_opt',	
					seed_name=	'wf1', mp_grid=	4, kubo_tol=1e-5,
					valence_bands=4,	
					n_hw=1,hw_min=0,hw_max=1.2,laser_phase=1,
					N_eF=301, eF_min=0, eF_max=3.0, Tkelvin=0,eta_smearing=0.1,
					plot_bands			=	False	, 
					debug_mode			=	False	,
                    use_cart_velo		=	False	, 
                    do_gauge_trafo		=	False	, 
                    R_vect_float		=	False	,
                    do_write_velo		=	False	, 
                    do_write_mep_bands	=	True	,
                    do_mep				=	True	, 
                    do_kubo				=	False	, 
                    do_ahc				=	False	, 
                    do_opt				=	False	, 
                    do_gyro				=	False
				):
		
		self.opt_log			=	[]
		#
		self.N_mpi_procs		=	N_mpi_procs
		#	DIRECTORY SETUP
		self.base_dir			= 	base_dir
		save_mkdir(base_dir)
		sh.copy('./mepInterp',	self.base_dir)
		#
		#	KSPACE
		self.seed_name			=	seed_name
		self.mp_grid			=	[mp_grid,mp_grid,mp_grid]
		self.kubo_tol			=	kubo_tol
		#
		#	MEP BANDS
		self.valence_bands		=	valence_bands
		#
		#	LASER
		self.n_hw				=	n_hw
		self.hw_min				=	hw_min
		self.hw_max				=	hw_max
		self.laser_phase		=	laser_phase
		#
		#	FERMI
		self.N_eF				=	N_eF
		self.eF_min				=	eF_min
		self.eF_max				=	eF_max
		self.Tkelvin			=	Tkelvin
		self.eta_smearing		=	eta_smearing
		#
		#	FLAGS
		self.plot_bands			=	plot_bands 
		self.debug_mode			=	debug_mode
		self.use_cart_velo		=	use_cart_velo 
		self.do_gauge_trafo		=	do_gauge_trafo 
		self.R_vect_float		=	R_vect_float
		self.do_write_velo		=	do_write_velo 
		self.do_write_mep_bands	=	do_write_mep_bands
		self.do_mep				=	do_mep 
		self.do_kubo			=	do_kubo 
		self.do_ahc				=	do_ahc 
		self.do_opt				=	do_opt 
		self.do_gyro			=	do_gyro


	def print_log(self):
		print("[MEP_optimizer]: History of the log	([t, strain, J_ex, tso,	theta_deg, 	cost]")
		for entry in self.opt_log:
			print(entry)

	def bandgap_exists(self,	occ_file, acc=1e-8):
		if os.path.isfile(occ_file):
			occ	=	np.genfromtxt(	occ_file, usecols=(1))
			print("len(occ)=",len(occ))
			for n_el in occ:
				#print("n_el=",n_el)
				if np.abs(	float(n_el - self.valence_bands)) < acc:
					return True
		else:
			print('[bandgap_exists]: ERROR occupation data file '+occ_file+' was not found!')
			return False
		return False

	def get_absmax_mep(self,	mep_file):
		mep_id		=	'mep'
		mep_tens	=	read_real_tens_file(	mep_file,	mep_id)
		mep_abs 	=	abs(	mep_tens	)
		
		print("mep_tens:",mep_tens)
		print("abs_max=",mep_abs.max())
		#
		return	mep_abs.max()







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def cost_func( x, args):	
		#	x:	t, strain, J_ex, tso, theta_deg
		print('[cost_func]: 	next x = ',x)
		costs		=	0
		#
		Jz			=	0
		#
		if x[0]	> 0: print('[cost_func]:	WARNING got corrupted para_val: t  =x[0]>0')
		if x[1]	< 0: print('[cost_func]:	WARNING got corrupted para_val: d  =x[1]<0')
		if x[2]	> 0: print('[cost_func]:	WARNING got corrupted para_val: J  =x[2]>0')
		if x[3]	> 0: print('[cost_func]:	WARNING got corrupted para_val: soc=x[3]>0')
		if x[4]	< 0: print('[cost_func]:	WARNING got corrupted para_val:  chi=x[4]<0')
		#
		new_run 	=	MEP_optimizer(	
							N_mpi_procs=4,
							mp_grid=32
						)
		#
		new_run_dir	=	'new_run'
		new_run_path=	new_run.base_dir+'/'+new_run_dir
		#save_mkdir(new_run.base_dir	+	//	new_run_dir)
		#
		write_3q_HR(	new_run.base_dir,new_run_dir, new_run.seed_name, 
						# optimizations paras:
						x[0],x[1], x[2],0,x[3], x[4],
						#	
	                    new_run.mp_grid, new_run.kubo_tol, new_run.valence_bands,
	                    new_run.n_hw, new_run.hw_min, new_run.hw_max,  new_run.laser_phase ,
	                    new_run.N_eF, new_run.eF_min, new_run.eF_max, new_run.Tkelvin,new_run.eta_smearing,
	                    new_run.plot_bands, new_run.debug_mode,
	                    new_run.use_cart_velo, new_run.do_gauge_trafo, new_run.R_vect_float,
	                    new_run.do_write_velo, new_run.do_write_mep_bands,
	                    new_run.do_mep, new_run.do_kubo, new_run.do_ahc, new_run.do_opt, new_run.do_gyro,
	                    False
	                )
		#
		#	run fortran code
		sh.copy(new_run.base_dir+'/mepInterp',	new_run_path+'/mepInterp')
		run_fortran(	new_run_path,	new_run.N_mpi_procs)
		#
		#
		occ_file	=	new_run_path 	+ 	'/out/occ.dat'
		mep_file	=	new_run_path	+	'/out/mep/mep_tens.dat'
		##	eval this run
		if(	new_run.bandgap_exists(occ_file)	):
			costs	=	-	new_run.get_absmax_mep(	mep_file	)
		else:
			print('[cost_func]: WARNING not an insulator')
			costs	=	0
		#
		print('{:16.8e}'.format(costs))
		return float(costs)

#def test():
#	new_run	=	MEP_optimizer(	mp_grid=16,	N_eF=301)
#	costs	=	new_run.cost_func(t=-1,strain=1.0,J_ex=-3.0,tso=0.0, theta_deg=	54.7)
#	#
#	new_run.print_log()
#test()




def run_opti(niter=100, verbose=False):
	init_guess			=	[-1.0,1.0,-3.0,-0.3,54.7]
	#
	args				=	[]
	constr				=	get_constraints()	
	#
	minimizer_kwargs	=	dict(method="COBYLA", constraints=constr, args=args)
	resHopp				=	basinhopping( cost_func,init_guess,	niter=niter, disp = verbose,	minimizer_kwargs=minimizer_kwargs)
	print("[run_opti]:	optimal parameter_set: ",resHopp.x)
	print("[run_opti]:	minimal value: ",resHopp.fun)
	print('[run_opti]:	optimizer terminated with final message:	')
	print('\t',resHopp.message)
#
run_opti(niter=1, verbose=True)
	






















