import os
import time
import datetime
import sys, traceback
import numpy as np
from shutil 			import 	copy
from tb_input_writer 	import 	write_souza_tb_input
from plot_bandStruct 	import	plot_bandstruct
from fortran_io			import	read_real_tens_file





class MEP_worker:
	"""class handeling a MEPinterp calculation"""


	#constructor
	def __init__(		self, root_dir, work_dir, phi, val_bands, mp_grid, 
						kubo_tol=1e-3,  hw=0.0, eFermi=0.0, Tkelvin=0.0, eta_smearing=0.0, 
						debug_mode='F', do_gauge_trafo='T',	
						do_write_velo='F', do_write_mep_bands='F',
						do_mep='T', do_kubo='F', do_ahc='F', do_opt='F', do_gyro='F'
						):
		self.root_dir		=	root_dir
		self.work_dir		=	work_dir
		self.phi			=	phi
		self.val_bands 		=	val_bands 
		self.mp_grid		=	mp_grid
		self.kubo_tol		=	kubo_tol
		self.hw				=	hw	
		self.eFermi			=	eFermi
		self.Tkelvin		=	Tkelvin
		self.eta_smearing	=	eta_smearing
		self.plot_bands		=	'F'	
		self.debug_mode		=	debug_mode
		self.do_gauge_trafo	=	do_gauge_trafo
		self.do_write_velo	=	do_write_velo
		self.do_write_mep_bands	=	do_write_mep_bands
		self.do_mep			=	do_mep
		self.do_kubo		=	do_kubo
		self.do_ahc			=	do_ahc
		self.do_opt			=	do_opt
		self.do_gyro		=	do_gyro
		#
		#
		self.success		= False
		self.band_dir		= self.work_dir+'/bands'
		#
		#
		print('************new calculation****************[   '+str(datetime.datetime.now())+'  ]')
		print('phi='+str(phi))
		print('mp_grid='+str(mp_grid))
		#make new folder
		try:
			os.makedirs(self.work_dir)
			print('made dir '+self.work_dir)
		except OSError:
			sys.exit('could not makedirs '+str(self.work_dir)	)

		#prepare files in working directory
		write_souza_tb_input(		self.work_dir, self.phi, self.val_bands, self.mp_grid, 					
									self.kubo_tol, self.hw, self.eFermi, self.Tkelvin,	self.eta_smearing,	 
									self.plot_bands, self.debug_mode, self.do_gauge_trafo,
									self.do_write_velo, self.do_write_mep_bands,
									self.do_mep, self.do_kubo, self.do_ahc, self.do_opt, self.do_gyro			 
							)
		copy(self.root_dir+'/../mepInterp',	self.work_dir)
		copy(self.root_dir+'/../kptsgen.pl',self.work_dir)





	#desctructor
	def __del__(self):
		print('done with calculation for phi='+str(self.phi) )
		print('--------------------------------------------------------')
		print('\n\n\n\n')
		try:
			os.rmdir(self.work_dir+'/raw')
		except:
			print('no "raw" dir found')
		finally:
			print('cleaned the working directory')




	#atributes
	def print_info(self):
		print('----current parameters:--------------')
		print('the root directory is',self.root_dir)
		print('the working directory is ',self.work_dir)
		print('phi= ',self.phi)
		print(' ')

	

	def run(self,mpi_np=1,plot_bands=False):
		#execute calculation
		
		os.chdir(self.work_dir)

		nK	= self.mp_grid[0]*self.mp_grid[1]*self.mp_grid[2]
		print('['+str(datetime.datetime.now())+']start calculation (nK='+str(nK)+')....')
		try:
			os.system('mpirun -np '+str(mpi_np)+' ./mepInterp > mep.log')
			self.success = True
		except:
			print('calculation could not be executed')

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
		return mep_tens, mep_cs, mep_lc, mep_ic

	def get_ahc_tens(self):
		ahc_file_path	= self.work_dir+'out/ahc/ahc_tens.dat'
		ahc_tens		= read_real_tens_file(ahc_file_path,					'ahc')
		return	ahc_tens

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
			print('could not make directory "'+self.band_dir+'" ')
		write_souza_tb_input(self.band_dir, self.phi, self.val_bands, self.mp_grid, 'T' )
		#
		#copy exectubales to target 
		copy(self.work_dir+'/mepInterp', 	self.band_dir)
		copy(self.work_dir+'/kptsgen.pl',	self.band_dir)
		#
		os.chdir(self.band_dir)
		print('['+str(datetime.datetime.now())+']start BAND calculation....')
		#genearte the kpt file
		os.system('./kptsgen.pl -l cub -k "Gamma 25 X 25 M 35 Gamma 25 R"')
		print('generated kpt file')
		#time.sleep(1)	#delay for one second

		#run the bandstructure calculation
		os.system('mpirun -np 4 ./mepInterp > mepBAND.log')
		print('MEPinterp run completed')

		#make a plot
		try:
			plot_bandstruct('kpts', 'MEPout/eBands.dat', 'bands.pdf',	label_size=14, y_tick_size=12, plot_in_ev=False)
		except:
			print('WARNING creation of bandstructure plot failed')

		
		try:
			os.rmdir(self.band_dir+'/raw')
		except OSError:
			print('could not remove "'+self.band_dir+'/raw", check for correct execution ')



		os.chdir(self.root_dir)
		print('['+datetime.datetime.now()+']...finished BAND plotting')

