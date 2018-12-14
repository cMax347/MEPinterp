import 	os
import sys
from shutil 			import 	copy
from shutil				import	rmtree
from tb_input_writer	import	write_tb_input
from plot_bandStruct	import	plot_bandstruct




def do_band_calc(tb_model, use_pos_op, latt_sym_group='cub',kpath="Gamma 1000 X 1000 M 1000 Gamma 1000 R", phi=0.0, val_bands=1, gauge_velos=True, write_velos=False):
	#create working directory
	root_dir	= os.getcwd()
	band_dir	= root_dir+'/bands'
	

	main_exe	= root_dir+'/mepInterp'
	kpt_gen		= root_dir+'/kptsgen.pl'

	
	#check if executables are present
	if os.path.isfile(main_exe) and os.path.isfile(kpt_gen):
		#delete old folder and create it again
		if os.path.isdir(band_dir):
			rmtree(band_dir)
			print('[do_band_calc]: [do_band_calc]:deleted old dir "',band_dir,'"')
		try:
			os.makedirs(band_dir)
		except OSError:
			print('[do_band_calc]: could not make directory "',band_dir,'"')
		#
		#now copy executables to target
		copy(main_exe,		band_dir+'/mepInterp')
		copy(kpt_gen,	band_dir+'/kptsgen.pl')
		print('[do_band_calc]: copied the executables')
		#
		#
		if tb_model == 'FeMn3q':
			inp_3q	= root_dir+'/inp_params_3q'
			if os.path.isfile(inp_3q):
				copy(inp_3q, band_dir+"/inp_params_3q")
				print("[do_band_calc]: copied '"+inp_3q+"' to "+band_dir)
			else:
				print("[do_band_calc]: could not find the 3q input file "+inp_3q)
				sys.exit()
	else:
		print('[do_band_calc]: did not find all executables necessary, nothing was done...')
		sys.exit()


	




	os.chdir(band_dir)

	#interpolate the energies
	mp_grid = [4,4,4]	#generic does not matter here

	#write_souza_tb_input(root_dir, phi_para, valence_bands, mp_grid , kubo_tol=1e-3, hw=0.0, eFermi=0.0, Tkelvin=0.0, eta_smearing=0.0, plot_bands='F'):

 

	write_tb_input(	tb_model, use_pos_op, band_dir, phi, val_bands, mp_grid, 
							kubo_tol=1e-3, 
							hw=0.0,laser_phase=0, 
							eFermi=0.0, Tkelvin=300.0, eta_smearing=3.0, 
							plot_bands='T', 
							debug_mode='T',
							do_gauge_trafo=gauge_velos, 
							R_vect_float='T',
							do_write_velo=write_velos
					)
	os.system('./kptsgen.pl -l '+latt_sym_group+' -k "'+kpath+'"')
	print('[do_band_calc]: generated k-space path list')
	os.system('mpirun -np 4 ./mepInterp > mepBANDS.log')
	print('[do_band_calc]: calculation done, now try plotting')


	#create pdf file 
	k_file	= band_dir+'/kpts'
	en_file	= band_dir+'/out/eBands.dat'
	pdf_file= band_dir+'/bands.pdf'
	if os.path.isfile(k_file) and os.path.isfile(en_file):
		plot_bandstruct(k_file,en_file, pdf_file, label_size=14, y_tick_size=12, plot_in_ev=True)
	else:
		print('[do_band_calc]: could not plot bandstructure since not all input files ("',k_file,'", "',en_file,'"") were found')


	os.chdir(root_dir)
	print('[do_band_calc]: by by')



do_band_calc(	tb_model		=	'FeMn3q'									,
				use_pos_op		=	False										,
				latt_sym_group	=	'jpH'										,
				kpath			=	'Gamma 1000 X 1000 M 1000 Gamma 1000 R'		, 
				phi				=	0.0											, 
				gauge_velos		=	False										,
				write_velos		=	True										,
			)


