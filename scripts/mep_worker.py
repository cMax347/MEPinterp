import os
import time
import datetime
import sys, traceback
import numpy as np
from shutil 			import 	copy
from tb_input_writer 	import 	write_tb_input
from tb_input_writer	import 	write_FeMn_3Q_inp
from plot_bandStruct 	import	plot_bandstruct
from fortran_io			import	read_real_tens_file





class MEP_worker:
	"""class handeling a MEPinterp calculation"""
	


	def write_3Q_input(self, t_hopp, delta, J_ex, lmbd_R):
		inp_file_3Q	=	self.work_dir+'/inp_params_3q'
		#
		write_FeMn_3Q_inp( inp_file_3Q	, t_hopp, delta, J_ex, lmbd_R)
		print("[MEP_worker]: wrote 3q input to "+inp_file_3Q)

	#constructor
	def __init__(		self, tb_model, use_pos_op, root_dir, work_dir, phi, val_bands, 
						t_hopp, delta, J_ex, lmbd_R,
						mp_grid, 
						kubo_tol=1e-3,  
						n_hw=1, hw_min=0.0, hw_max=1.0, laser_phase=1, 
						n_eF=1 ,eF_min=0.0, eF_max=0.0, Tkelvin=0.0, eta_smearing=0.0, 
						debug_mode='F', do_gauge_trafo='T',	R_vect_float='F',
						do_write_velo='F', do_write_mep_bands='F',
						do_mep='T', do_kubo='F', do_ahc='F', do_opt='F', do_gyro='F'
						):
		self.tb_model		=	tb_model
		self.use_pos_op		=	use_pos_op
		self.root_dir		=	root_dir
		self.work_dir		=	work_dir
		#
		self.phi			=	phi
		self.val_bands 		=	val_bands
		#
		self.t_hopp  		=	t_hopp	
		self.delta   		=	delta
		self.J_ex    		=	J_ex
		self.lmbd_R  		=	lmbd_R
		#
		self.mp_grid		=	mp_grid
		self.kubo_tol		=	kubo_tol
		#
		self.n_hw			=	n_hw	
		self.hw_min			=	hw_min
		self.hw_max			=	hw_max		
		#
		self.laser_phase	=	laser_phase
		self.n_eF			=	n_eF
		self.eF_min			=	eF_min
		self.eF_max			=	eF_max
		self.Tkelvin		=	Tkelvin
		self.eta_smearing	=	eta_smearing
		#
		self.plot_bands		=	'F'	
		self.debug_mode		=	debug_mode
		self.do_gauge_trafo	=	do_gauge_trafo
		self.R_vect_float	=	R_vect_float
		self.do_write_velo	=	do_write_velo
		self.do_write_mep_bands	=	do_write_mep_bands
		#
		self.do_mep			=	do_mep
		self.do_kubo		=	do_kubo
		self.do_ahc			=	do_ahc
		self.do_opt			=	do_opt
		self.do_gyro		=	do_gyro
		#
		#
		self.success		= False
		print('[MEP_worker/init]: work dir=',self.work_dir)
		self.band_dir		= self.work_dir+'/bands'
		#
		#
		print('[mep_worker]: ************new model system****************[   '+str(datetime.datetime.now())+'  ]')
		#print('phi='+str(phi))
		print('[mep_worker]: mp_grid='+str(mp_grid))
		#make new folder
		try:
			os.makedirs(self.work_dir)
			print('[mep_worker]: made dir '+self.work_dir)
		except OSError:
			print('[mep_worker]: could not makedirs '+str(self.work_dir)	)

		#copy executables to target
		copy(self.root_dir+'/../mepInterp',	self.work_dir)
		copy(self.root_dir+'/../kptsgen.pl',self.work_dir)
		
		self.write_3Q_input(	t_hopp, delta, J_ex, lmbd_R	)
		## copy existing input
		#if self.tb_model=='FeMn3q':
		#	#try:
		#	#	copy(self.root_dir+'/../inp_params_3q',	self.work_dir)
		#	#except:
		#	#	print('ERROR file '+self.root_dir+'/../inp_params_3q not found')

		#generate fortran input
		write_tb_input(		self.tb_model, self.use_pos_op,	self.work_dir, 
							self.phi, self.val_bands, self.mp_grid, 					
							self.kubo_tol, 
							self.n_hw, self.hw_min, self.hw_max , self.laser_phase, 
							self.n_eF, self.eF_min, self.eF_max, self.Tkelvin,	self.eta_smearing,	 
							self.plot_bands, self.debug_mode, self.do_gauge_trafo, self.R_vect_float,
							self.do_write_velo, self.do_write_mep_bands,
							self.do_mep, self.do_kubo, self.do_ahc, self.do_opt, self.do_gyro			 
							)





	

	#atributes
	def print_info(self):
		print('[mep_worker]: ----infos about current run:--------------')
		print('[mep_worker]: the root directory is',self.root_dir)
		print('[mep_worker]: the working directory is ',self.work_dir)
		#print('phi= ',self.phi)
		print(' ')

	

	def run(self, dry=False, mpi_np=1,plot_bands=False):
		#execute calculation
		
		os.chdir(self.work_dir)

		nK	= self.mp_grid[0]*self.mp_grid[1]*self.mp_grid[2]
		print('['+str(datetime.datetime.now())+']start calculation (nK='+str(nK)+')....')
		if not dry:
			try:
				os.system('mpirun -np '+str(mpi_np)+' ./mepInterp > mepINTERPOLATION.log')
				self.success = True
			except:
				print('[mep_worker]: calculation could not be executed')
		else:
			print('[mep_worker]: dry run selected no computation done!')

		os.chdir(self.root_dir)
		print('['+str(datetime.datetime.now())+'] ...finished calculation')
		


	def get_mep_tens(self):
		#TOTAL
		mep_file_path	= self.work_dir+'/out/mep/mep_tens.dat'
		mep_tens		= read_real_tens_file(mep_file_path,					'mep')
		#CHERN-SIMONS
		mep_file_path	= self.work_dir+'/out/mep/mep_cs.dat'
		mep_cs			= read_real_tens_file(mep_file_path,					'mep')
		#LOCAL
		mep_file_path	= self.work_dir+'/out/mep/mep_lc.dat'
		mep_lc			= read_real_tens_file(mep_file_path,					'mep')
		#ITINERANT
		mep_file_path	= self.work_dir+'/out/mep/mep_ic.dat'
		mep_ic			= read_real_tens_file(mep_file_path,					'mep')
		#
		#
		#BAND RESOLVED
		mep_bands 		= []
		for n in range(1,self.val_bands+1):
			mep_file_path	= self.work_dir+'/out/mep/mep_band.'+"{:07d}".format(n)
			tmp	=	read_real_tens_file(mep_file_path,				'mep')
			mep_bands.append(tmp)



		return mep_tens, mep_cs, mep_lc, mep_ic, mep_bands

	def get_ahc_tens(self):
		ahc_file_path	= self.work_dir+'out/ahc/ahc_tens.dat'
		ahc_tens		= read_real_tens_file(ahc_file_path,					'ahc')
		return	ahc_tens


	def get_ohc_tens(self):
		ohc_file_path	=	self.work_dir+'out/ahc/ahc_tens.dat'

	def get_opt_tens(self):
		#	symmetric contribution
		Ssymm_file_path	=	self.work_dir+'out/opt/opt_Ssymm.dat'
		Ssymm_tens		=	read_real_tens_file(Ssymm_file_path,				'optS')
		#
		#	a symmetric contribution
		Asymm_file_path	=	self.work_dir+'out/opt/opt_Asymm.dat'
		Asymm_tens		=	read_real_tens_file(Asymm_file_path,				'optA')
		#
		#
		return Ssymm_tens, Asymm_tens


	def plot_bands(self):
		#
		try:
			os.mkdir(self.band_dir)
		except OSError:
			print('[mep_worker]: could not make directory "'+self.band_dir+'" ')
		write_souza_tb_input(self.band_dir, self.phi, self.val_bands, self.mp_grid, 'T' )
		#
		#copy exectubales to target 
		copy(self.work_dir+'/mepInterp', 	self.band_dir)
		copy(self.work_dir+'/kptsgen.pl',	self.band_dir)

		#
		if self.tb_model=="FeMn3q":
			copy(self.work_dir+'/inp_params_3q',	self.band_dir)
		#
		os.chdir(self.band_dir)
		print('['+str(datetime.datetime.now())+';mep_worker]: start BAND calculation....')
		#genearte the kpt file
		os.system('./kptsgen.pl -l cub -k "Gamma 25 X 25 M 35 Gamma 25 R"')
		print('generated kpt file')
		#time.sleep(1)	#delay for one second

		#run the bandstructure calculation
		os.system('mpirun -np 4 ./mepInterp > mepBAND.log 2> mepBAND.err')
		print('[mep_worker]: MEPinterp run completed')

		#make a plot
		try:
			plot_bandstruct('kpts', 'MEPout/eBands.dat', 'bands.pdf',	label_size=14, y_tick_size=12, plot_in_ev=False)
		except:
			print('[mep_worker]: WARNING creation of bandstructure plot failed')

		
		try:
			os.rmdir(self.band_dir+'/raw')
		except OSError:
			print('[mep_worker]: could not remove "'+self.band_dir+'/raw", check for correct execution ')



		os.chdir(self.root_dir)
		print('['+datetime.datetime.now()+';mep_worker]: ...finished BAND plotting')

