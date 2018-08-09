import 	os
import sys
from shutil 			import 	copy
from shutil				import	rmtree
from tb_input_writer	import	write_souza_tb_input
from plot_bandStruct	import	plot_bandstruct




def do_band_calc(phi, val_bands):
	#create working directory
	root_dir	= os.getcwd()
	band_dir	= root_dir+'/bands'
	

	main_exe	= root_dir+'/main.exe'
	kpt_gen		= root_dir+'/kptsgen.pl'

	#check if executables are present
	if os.path.isfile(main_exe) and os.path.isfile(kpt_gen):
		#delete old folder and create it again
		if os.path.isdir(band_dir):
			rmtree(band_dir)
			print('deleted old dir "',band_dir,'"')
		try:
			os.makedirs(band_dir)
		except OSError:
			print('could not make directory "',band_dir,'"')
		#
		#now copy executables to target
		copy('./main.exe',		band_dir+'/main.exe')
		copy('./kptsgen.pl',	band_dir)
		print('copied the executables')
	else:
		print('did not find all executables necessary, nothing was done...')
		sys.exit()

	




	os.chdir(band_dir)

	#interpolate the energies
	mp_grid = [4,4,4]	#generic does not matter here
	write_souza_tb_input(band_dir, phi, val_bands, mp_grid, use_interp_kpt='F', do_gauge_trafo='T', plot_bands='T' )
	os.system('./kptsgen.pl -l cub -k "Gamma 25 X 25 M 35 Gamma 35 R"')
	print('generated k-space path list')
	os.system('mpirun -np 4 ./main.exe > mepBANDS.out')
	print('calculation done, now try plotting')


	#create pdf file 
	k_file	= band_dir+'/kpts'
	en_file	= band_dir+'/MEPout/eBands.dat'
	pdf_file= band_dir+'/bands.pdf'
	if os.path.isfile(k_file) and os.path.isfile(en_file):
		plot_bandstruct(k_file,en_file, pdf_file, label_size=14, y_tick_size=12, plot_in_ev=False)
	else:
		print('could not plot bandstructure since not all input files ("',k_file,'", "',en_file,'"") were found')


	os.chdir(root_dir)
	print('by by')



do_band_calc(0.0,2)

