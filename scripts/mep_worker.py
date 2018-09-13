import os
import time
import datetime
import sys, traceback
import numpy as np
from shutil 			import 	copy
from tb_input_writer 	import 	write_souza_tb_input
from plot_bandStruct 	import	plot_bandstruct
from fortran_io			import 	read_mep_file





class MEP_worker:
	"""class handeling a MEPinterp calculation"""


	#constructor
	def __init__(self, root_dir, work_dir, phi, val_bands, mp_grid, use_interp_kpt='F', do_gauge_trafo='T'	):
		self.root_dir		= root_dir
		self.work_dir		= work_dir
		self.phi			= phi
		self.val_bands 		= val_bands 
		self.mp_grid		= mp_grid
		self.use_interp_kpt = use_interp_kpt
		self.do_gauge_trafo	= do_gauge_trafo
		#
		self.success		= False
		self.band_dir		= self.work_dir+'/bands'
		#
		#
		print('************new calculation****************')
		print('phi='+str(phi))
		print('mp_grid='+str(mp_grid))
		#make new folder
		try:
			os.makedirs(self.work_dir)
			print('made dir '+self.work_dir)
		except OSError:
			sys.exit('could not makedirs '+str(self.work_dir)	)

		#prepare files in working directory
		write_souza_tb_input(self.work_dir, self.phi, self.val_bands, self.mp_grid, self.use_interp_kpt, self.do_gauge_trafo, 'F' )
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
		mep_file_path	= self.work_dir+'/mep/mep_tens.dat'
		mep_tens		= read_mep_file(mep_file_path)
		#CHERN-SIMONS
		mep_file_path	= self.work_dir+'/mep/mep_cs.dat'
		mep_cs			= read_mep_file(mep_file_path)
		#LOCAL
		mep_file_path	= self.work_dir+'/mep/mep_lc.dat'
		mep_lc		= read_mep_file(mep_file_path)
		#ITINERANT
		mep_file_path	= self.work_dir+'/mep/mep_ic.dat'
		mep_ic		= read_mep_file(mep_file_path)


		return mep_tens, mep_cs, mep_lc, mep_ic




	def plot_bands(self):
		#
		try:
			os.mkdir(self.band_dir)
		except OSError:
			print('could not make directory "'+self.band_dir+'" ')
		write_souza_tb_input(self.band_dir, self.phi, self.val_bands, self.mp_grid, self.use_interp_kpt, self.do_gauge_trafo, 'T' )
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

