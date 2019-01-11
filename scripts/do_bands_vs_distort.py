import 	os
import sys
import datetime
import numpy 			as		np
from shutil 			import 	copy
from shutil				import	rmtree
from tb_input_writer	import	write_tb_input
from tb_input_writer	import	write_FeMn_3Q_inp
from plot_bandStruct	import	plot_bandstruct



def save_mkdir(dir):
	if os.path.isdir(dir):
		rmtree(dir)
		print('[save_mkdir]: [do_band_calc]:deleted old dir "',dir,'"')
	try:
		os.makedirs(dir)
		print('[ave_mkdir]: I created '+dir)
	except OSError:
		print('[save_mkdir]: could not make directory "',dir,'"')



def save_root_dir_creator(orig_root, root_name):
	#derived attributes
	root_dir	= orig_root+'/'+root_name+'_'+datetime.date.today().strftime("%d%B%Y")
	#search new root folder name
	#create root directory
	if os.path.isdir(root_dir):
		old_path 	= root_dir
		cnt 		= 0
		while os.path.isdir(old_path) and cnt < 5:
			old_path	= old_path + '.old'
			cnt			= cnt + 1
			try:
				os.rename(root_dir, old_path)
			except OSError:
				print(old_path+ ' exists already')
	if os.path.isdir(root_dir):
		stop
	os.mkdir(root_dir)
	#
	#
	return root_dir



def prepare_calc(		target_dir, tb_model, use_pos_op, val_bands, main_exe, kpt_gen, J_ex, t1_hopp, delta,
						gauge_velos, write_velos
				):
	copy(	main_exe	,		target_dir	+	'/mepInterp')
	copy(	kpt_gen		,		target_dir	+	'/kptsgen.pl')
	print('[prepare_calc]: copied the executables')
	#	write 3q input file
	if tb_model == 'FeMn3q':
		#	
		#	write input file
		inp_3q	= target_dir	+	'/inp_params_3q'
		write_FeMn_3Q_inp(inp_3q, t1_hopp, delta, J_ex)
		#
		#interpolate the energies
		mp_grid = [4,4,4]	#generic does not matter here
		phi = 0.0			#generic does not matter here
		write_tb_input(	tb_model, use_pos_op, target_dir, phi, val_bands, mp_grid, 
								kubo_tol=1e-3, 
								hw=0.0,laser_phase=0, 
								 N_eF=1, eF_min=0.0, eF_max=0, Tkelvin=300.0, eta_smearing=3.0, 
								plot_bands='T', 
								debug_mode='T',
								do_gauge_trafo=gauge_velos, 
								R_vect_float='T',
								do_write_velo=write_velos
						)


def run_job(root_dir, job_dir, latt_sym_group, kpath):
	os.chdir(job_dir)
	#
	os.system('./kptsgen.pl -l '+latt_sym_group+' -k "'+kpath+'"')
	print('[run_job]: generated k-space path list')
	os.system('mpirun -np 4 ./mepInterp > mepBANDS.log')
	k_file	= job_dir+'/kpts'
	en_file	= job_dir+'/out/eBands.dat'
	#
	os.chdir(root_dir)
	print('[run_job]: calculation done')
	#
	return k_file, en_file





def sample_distortion(	tb_model, use_pos_op, 
						latt_sym_group='cub',
						kpath="Gamma 1000 X 1000 M 1000 Gamma 1000 R", 
						val_bands=1, gauge_velos=True, write_velos=False,
						J_ex			=	1.0											,
						t_hopp			=	1.0											,
						delta_min	=	0.9,	delta_max	=	1.0,
						delta_steps	=	5
					):
	

	print('[sample_distortion]: hello there')
	#create working directory

	
	orig_root	=	os.getcwd()
	root_dir	=	save_root_dir_creator(orig_root,'delta_band_run')
	main_exe	= 	root_dir+'/mepInterp'
	kpt_gen		= 	root_dir+'/kptsgen.pl'
	#
	copy(	orig_root	+	'/mepInterp'	,	main_exe	)
	copy(	orig_root	+	'/kptsgen.pl'	, 	kpt_gen		)
	#
	job_dirs	=	[]
	delta_lst	=	[]
	k_file_lst	=	[]
	en_file_lst	=	[]
	#
	#
	for delta in np.linspace(delta_min,	delta_max, delta_steps):
		delta	=	round(delta,3)
		print("[sample_distortion]:	start job delta="+str(delta))
		band_dir	= root_dir+'/bands_d'+str(delta)
		job_dirs.append(	band_dir )
		#
		#check if executables are present
		if os.path.isfile(main_exe) and os.path.isfile(kpt_gen):
			#delete old folder and create it again
			save_mkdir(band_dir)
			#
			#now copy executables to target
			prepare_calc(band_dir, tb_model, use_pos_op, val_bands, main_exe, kpt_gen, J_ex, t_hopp, delta, gauge_velos, write_velos)
			#
			#run exe in target
			new_k, new_en	=	run_job(root_dir, band_dir, latt_sym_group, kpath)
			#
			delta_lst.append(	delta	)
			k_file_lst.append(	new_k	)
			en_file_lst.append(	new_en	)
			#
		else:
			print('[sample_distortion]: did not find all executables necessary, nothing was done...')
			sys.exit()
		#
		#
	return delta_lst, k_file_lst, en_file_lst, job_dirs

	



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	




	##create pdf file 
#	#pdf_file= band_dir+'/bands.pdf'
#	#if os.path.isfile(k_file) and os.path.isfile(en_file):
#	#	plot_bandstruct(k_file,en_file, pdf_file, label_size=14, y_tick_size=12, plot_in_ev=True)
#	#else:
#	#	print('[do_band_calc]: could not plot bandstructure since not all input files ("',k_file,'", "',en_file,'"") were found')
#
#
	#os.chdir(root_dir)
	#print('[do_band_calc]: by by')



def main():
	pdf_file	=	'./distorted_bands.pdf'
	id_str		=	'bands_d'
	id_formula	=	r'$\delta$'
	#
	delta_lst, k_file_lst, en_file_lst,	job_dirs	=	sample_distortion(	
					tb_model		=	'FeMn3q'									,
					use_pos_op		=	False										,
					latt_sym_group	=	'jpH'										,
					kpath			=	'Gamma 1000 X 1000 M 1000 Gamma 1000 R'		, 
					gauge_velos		=	False										,
					write_velos		=	False										,
					J_ex			=	1.0											,
					t_hopp			=	1.0											,
					delta_min		=	0.9											,
					delta_max		=	1.1											,
					delta_steps		=	5
				)
	print("the following job dirs where returned: ",job_dirs)
	print("\n")
	print("returned energy file list:")
	for en_file in en_file_lst:
		print(en_file)
	print("returned delta list:")
	for delta in delta_lst:
		print(delta)

	print("now try to plot everyhting \n ...")


	line_style	=	['-',':','-',':','-']
	plot_color	=	['orange','orange','black','blue','blue']


	plot_bandstruct(	job_dirs, id_str, id_formula, line_style, plot_color, pdf_file, label_size=14, y_tick_size=12, plot_in_ev=True )
#
main()


