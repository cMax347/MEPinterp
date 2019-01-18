import 	os
import 	datetime
import 	numpy 				as 		np
from 	shutil 				import 	rmtree
from 	postw90_in_writer 	import 	postw90_job
from 	souza_tb_model		import 	get_souza_tb
from 	FeMn_3q_model		import 	get_FeMn3q_tb
#
#
#
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#		CONSTANTS
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
whitespace 		= 	'\t'
#
au_to_eV 		= 	27.21139
au_to_angtrom	= 	0.529177211
#
seed_name		= 	'wf1'
#
#
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#		WANNIER90 INTERFACE
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def write_postw90_input(w90_dir, seed_name, val_bands, mp_grid, hw=0.0, eFermi=0.0, Tkelvin=0.0, eta_smearing=0.0):
	pw90_job 	=	postw90_job(w90_dir, seed_name, val_bands, mp_grid,  hw, eFermi, Tkelvin, eta_smearing)
	pw90_job.write_win_file()
#
#
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#		WRITE TIGHT BINDING BASIS (_hr, _r files)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def write_hr_file(seed_name, nAt, nrpts, thopp):
	
	#
	#print("thopp")
	#for elem in thopp:
	#	print(elem)
	thopp = sorted(thopp)
	with open(seed_name+'_hr.dat','w') as outfile:
		#HEADER
		outfile.write('# my _hr file created with inpgenTBwann.py on '+datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")+'\n')
		outfile.write(str(nAt)+'\n')
		outfile.write(str(nrpts)+'\n')
		#WIGNER SEITZ DEGNERACY
		cnt = 1
		for rpt in range(nrpts):
			outfile.write(str(2)+' ')
			cnt = cnt +1
			if cnt is 15:
				outfile.write('\n')
				cnt = 1
		outfile.write('\n')
		#
		#BODY
		for t in thopp:
			#write integers
			for val in t[0:5]:
				outfile.write(str(val)+whitespace)
			#write floats
			for val in t[5:7]:
				outfile.write("{:16.8f}".format(float(val))+whitespace)
			outfile.write('\n')
#
#
def write_r_file(seed_name, nAt, rhopp ):
	whitespace = '\t'
	rhopp = sorted(rhopp)
	#
	with open(seed_name+'_r.dat','w') as outfile:
		#HEADER
		outfile.write('# my _r file created with inpgenTBwann.py on '+datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")+'\n')
		outfile.write(str(nAt)+'\n')
		#
		#BODY
		for r in rhopp:
			#integers
			for val in r[0:5]:
				outfile.write(str(val)+whitespace)
			#floats
			for val in r[5:11]:
				outfile.write("{:16.8f}".format(float(val))+whitespace)
			outfile.write('\n')

#
#
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#		CFG FILE	WRITERS
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def write_mepInterp_input(	file_path,valence_bands, ax, ay, az, a0, mp_grid, seed_name,			
							kubo_tol, n_hw, hw_min, hw_max,  laser_phase ,N_eF, eF_min, eF_max, Tkelvin,eta_smearing,  	 
							plot_bands,	debug_mode, do_gauge_trafo, R_vect_float	, do_write_velo,	do_write_mep_bands,								
							do_mep, do_kubo, do_ahc, do_opt, do_gyro		
						):
	with open(file_path+'input.cfg','w') as outfile:
		outfile.write('# input file for TB model from New J Physics 12, 053032 (2010)'+'\n')
		outfile.write('# generated on '+datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")+'\n')
		outfile.write('\n')
		#
		#
		outfile.write('[jobs]\n')
		outfile.write('    '	+	'plot_bands='		+	str(plot_bands)			+	'\n')
		outfile.write('    '	+	'debug_mode='		+	str(debug_mode)			+	'\n')		
		outfile.write('    '	+	'R_vect_float='		+	str(R_vect_float)		+	'\n')
		outfile.write('    '	+	'do_write_velo='	+	str(do_write_velo)		+	'\n')		
		outfile.write('    '	+	'do_mep='			+	str(do_mep)				+	'\n')		
		outfile.write('    '	+	'do_kubo='			+	str(do_kubo)			+	'\n')		
		outfile.write('    '	+	'do_ahc='			+	str(do_ahc)				+	'\n')		
		outfile.write('    '	+	'do_opt='			+	str(do_opt)				+	'\n')		
		outfile.write('    '	+	'do_gyro='			+	str(do_gyro)			+	'\n')		
		outfile.write('\n')
		#
		#
		outfile.write('[unitCell]\n')
		outfile.write('    '	+	'a1= '				+	str(ax[0]) +' '+str(ax[1])+' '+str(ax[2])		+	'\n')
		outfile.write('    '	+	'a2= '				+	str(ay[0]) +' '+str(ay[1])+' '+str(ay[2])		+	'\n')
		outfile.write('    '	+	'a3= '				+	str(az[0]) +' '+str(az[1])+' '+str(az[2])		+	'\n')
		outfile.write('    '	+	'a0= '				+	str(a0)				+	'\n')
		outfile.write('\n')
		#
		#
		outfile.write('[wannInterp]\n')
		outfile.write('    '	+	'doGaugeTrafo= '	+	str(do_gauge_trafo)		+	'\n')
		outfile.write('    '	+	'mp_grid= '			+	str(mp_grid[0]) +' '+str(mp_grid[1])+' '+str(mp_grid[2])		+	'\n')
		outfile.write('    '	+	'seed_name= '		+	seed_name			+	'\n')
		outfile.write('\n')
		#
		#
		outfile.write('[MEP]\n')
		outfile.write('    '	+	'valence_bands= '	+	str(valence_bands)	+	'\n')
		outfile.write('    '	+	'do_write_mep_bands= '+	str(do_write_mep_bands)+'\n')
		outfile.write('\n')
		#
		#
		outfile.write('[Fermi]\n')
		outfile.write('    '	+	'N_eF= '			+	str(N_eF)			+	'\n')
		outfile.write('    '	+	'eF_min= '			+	str(eF_min)			+	'\n')
		outfile.write('    '	+	'eF_max= '			+	str(eF_max)			+	'\n')
		outfile.write('    '	+	'Tkelvin= '			+	str(Tkelvin)		+	'\n')
		outfile.write('    '	+	'eta_smearing= '	+	str(eta_smearing)	+	'\n')
		outfile.write('    '	+	'kuboTol= '			+	str(kubo_tol)		+	'\n')
		outfile.write('\n')
		#
		#


		#call CFG_add_get(my_cfg,	"Laser%N_hw"					,	N_hw					,	"points to probe in interval"	)
		#		call CFG_add_get(my_cfg,	"Laser%hw_min"					,	hw_min					,	"min energy of incoming light"	)
		#		call CFG_add_get(my_cfg,	"Laser%hw_max"					,	hw_max					,	"max energy of incoming light"	)
		#		call CFG_add_get(my_cfg,	"Laser%laser_phase"				,	laser_phase			,	"euler angle of phase shift of driving E-field")
		#		!
		outfile.write('[Laser]\n')
		outfile.write('    '	+	'N_hw= '			+	str(int(n_hw))  	+	'\n')
		outfile.write('    '	+	'hw_min= '			+	str(hw_min)  		+	'\n')
		outfile.write('    '	+	'hw_max= '			+	str(hw_max)  		+	'\n')
		outfile.write('    '	+	'laser_phase= '		+	str(laser_phase)	+	'\n')



def write_FeMn_3Q_inp(
				inp3Q_file, t1_hopp, delta, J_ex, lmbd_R=0,
				phiA	=	0.0		,		thetaA	=	0.0,
				phiB	=	30.0	,		thetaB	=	-109.47122063449069,
				phiC	=	270.0	,		thetaC	=	-109.47122063449069,
				phiD	=	150.0	,		thetaD	=	-109.47122063449069
					):
	#
	if os.path.isfile(inp3Q_file):	os.remove(inp3Q_file)
	#
	#	setup interlayer hopping t2_hopp
	if np.abs( (delta)**2 ) > 1e-4:
		t2_hopp	= t1_hopp / (delta)**2
	else:
		t2_hopp	= t1_hopp
	#
	hopp_str	=	str(t1_hopp)+" "+str(t2_hopp)+" "+str(J_ex)+"		"+"    ! t1, t2, lambda\n"
	#
	#	write inp_params_3q file
	with open(inp3Q_file,'w') as config_file:
		config_file.write("16 16 16                     ! Nkx, Nky, Nkz for BZ integration	\n")
		config_file.write(hopp_str)
		config_file.write(	"   "	+	str(phiA)	+	"   "	+	str(thetaA) +	"               ! phiA, thetaA for spin A	\n")
		config_file.write(	"   "	+	str(phiB)	+	"   "	+	str(thetaB) +	"               ! phiB, thetaB for spin B	\n")
		config_file.write(	"   "	+	str(phiC)	+	"   "	+	str(thetaC) +	"               ! phiC, thetaC for spin C	\n")
		config_file.write(	"   "	+	str(phiD)	+	"   "	+	str(thetaD) +	"               ! phiD, thetaD for spin A	\n")
		config_file.write(" -5.0 -3.0  200               ! Efmin, Efmax, NEf		\n")
		config_file.write(" F                            ! k-path instead?			\n")
		config_file.write(" 2000                         ! Nk for k-path			\n")
		config_file.write(	"   "	+	str(lmbd_R)+	"               ! Rashba SOC strength in eV (additional para by max merte)\n")

#
#
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#		PUBLIC FUNCTION
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def write_tb_input(	tb_model,use_pos_op, root_dir, phi_para, valence_bands, mp_grid ,						
 							kubo_tol, 
 							n_hw, hw_min, hw_max, laser_phase,  
 							N_eF, eF_min, eF_max, Tkelvin, eta_smearing, 	
 							plot_bands='F', debug_mode='F' ,do_gauge_trafo='T',	 R_vect_float='F',
 							do_write_velo='F',	do_write_mep_bands='F',			
 							do_mep='T', do_kubo='F', do_ahc='F', do_opt='F', do_gyro='F'		
 						):
	target_dir_name	= 'w90files'
	target_path		= root_dir+'/'+target_dir_name
	#
	#		get the tight binding basis
	if tb_model == 'souza':
		print("[write_tb_input]:	SOUZA TB MODEL SELECTED")
		nWfs, nrpts, tHopp, rHopp, ax, ay, az, a0	=	get_souza_tb(phi_para)
	elif tb_model == 'FeMn3q':
		print("[write_tb_input]:	FeMn 3Q TB MODEL SELECTED")
		config_file =	root_dir+'/inp_params_3q'
		nWfs, nrpts, tHopp, rHopp, ax, ay, az, a0	=	get_FeMn3q_tb(config_file, verbose=True)
	elif tb_model == 'w90':
			print("[write_tb_input]:	USE EXISTING W90 FILES (NOT IMPLEMENTED YET)")
			#	todo check if w90 files are present & consistent
	else:
		print("[write_tb_input]:	unknow tb_model identifier '"+tb_model+"' will use souza tb model with phi=0.0")
		nWfs, nrpts, tHopp, rHopp, ax, ay, az, a0	=	get_souza_tb(0.0)
	#
	#
	#write them to file
	os.mkdir(target_path)
	write_hr_file(	target_path+'/'+seed_name,	nWfs, nrpts, tHopp )
	print("[write_tb_input]:	wrote "+target_path+'/'+seed_name+'_hr.dat')
	if use_pos_op:
		write_r_file(	target_path+'/'+seed_name, 	nWfs,		 rHopp )
		print("[write_tb_input]:	wrote "+target_path+'/'+seed_name+'_r.dat')

	# now write the input files for postw90 & mepInterp
	#write_postw90_input(target_path, seed_name, valence_bands, mp_grid,  hw, eFermi, Tkelvin, eta_smearing	)
	#print("[write_tb_input]: wrote pw90 input file")

	#
	
	#	write mepInterp
	write_mepInterp_input(	 root_dir+'/',valence_bands, ax, ay, az, a0, mp_grid, seed_name,	
							kubo_tol, n_hw, hw_min, hw_max,  laser_phase, N_eF, eF_min, eF_max, Tkelvin,eta_smearing, 							
							plot_bands, debug_mode, do_gauge_trafo,	R_vect_float,
							do_write_velo, do_write_mep_bands,							
							do_mep, do_kubo, do_ahc, do_opt, do_gyro							
						)
	print("[write_tb_input]: wrote "+str(root_dir)+"/input.cfg file for fortran core")

#
#
#def write_FeMn3q_tb_input():
#	print("[write_FeMn3q_tb_input]:	Implement me!!!")




#
#
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#		TEST FUNCTION FOR DEBUGING ETC.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def test():
	print('+++++++++++++++++++++++++INPGEN-TB-WANN+++++++++++++++++++++++++++')
	print('generates wannier90 style input files based on the TB model in New J Physics 12, 053032 (2010)')
	root_dir = "test_tb_input_writer"
	print('test the tb input generator, files will be written within ',root_dir)
	if os.path.isdir(root_dir):
		rmtree(root_dir)
		print('removed old dir ',root_dir)
	os.mkdir(root_dir)
	print('created diretory ',root_dir)

	phi_para = 0.0
	valence_bands=2
	mp_grid=[4, 4, 4]

	write_tb_input('souza', root_dir, phi_para, valence_bands, mp_grid)


#uncomment next line to generate some test input
#test()


