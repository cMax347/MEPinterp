import os
import time
import datetime
import sys, traceback
import numpy 			as 		np
from shutil 			import 	copy
from tb_input_writer 	import 	write_souza_tb_input
from plot_bandStruct 	import	plot_bandstruct





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
		copy(self.root_dir+'/../main.exe',	self.work_dir)
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

	

	def run(self,plot_bands=False):
		#execute calculation
		
		os.chdir(self.work_dir)

		print('['+str(datetime.datetime.now())+']start calculation....')
		try:
			os.system('mpirun -np 4 ./main.exe > mep.log')
			self.success = True
		except:
			print('calculation could not be executed')

		os.chdir(self.root_dir)
		print('['+str(datetime.datetime.now())+'] ...finished calculation')
		


	def get_mep_tens(self):
		mep_tens = []
		if self.success:
			mep_file_path	= self.work_dir+'/MEPout/mep_tens.dat'
			with open(mep_file_path, 'r') as mep_file:
				start = -10
				for idx,line in enumerate(mep_file):
					if 'begin mep' in line:
						start = idx
					if idx > start  and idx <= start + 3:
						mep_tens.append(np.fromstring( line, dtype=np.float, sep=' ' ) )
					if idx == start + 4 and 'end mep' not in line:
						print('error at the end of reading mep tensor')

			mep_tens = np.array(mep_tens)
			print('read file '+mep_file_path)
			if mep_tens.size is not 9:
				print('WARNING issues with dimensionalty of MEP tensor')
				print('mep_tens interpretation: '+str(mep_tens))
		else:
			print('calculation was not successfull, can not grep MEP tensor.')

		return mep_tens





	def plot_bands(self):
		#
		try:
			os.mkdir(self.band_dir)
		except OSError:
			print('could not make directory "'+self.band_dir+'" ')
		write_souza_tb_input(self.band_dir, self.phi, self.val_bands, self.mp_grid, self.use_interp_kpt, self.do_gauge_trafo, 'T' )
		#
		#copy exectubales to target 
		copy(self.work_dir+'/main.exe', 	self.band_dir)
		copy(self.work_dir+'/kptsgen.pl',	self.band_dir)
		#
		os.chdir(self.band_dir)
		print('['+str(datetime.datetime.now())+']start BAND calculation....')
		#genearte the kpt file
		os.system('./kptsgen.pl -l cub -k Gamma 25 X 25 M 35 Gamma 25 R')
		print('generated kpt file')
		#time.sleep(1)	#delay for one second

		#run the bandstructure calculation
		os.system('mpirun -np 4 ./main.exe > mepBAND.log')
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

