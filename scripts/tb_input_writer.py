import numpy as np
import datetime
import os
from	shutil import rmtree
from postw90_in_writer import postw90_job
from souza_tb_model	import get_souza_tb


au_to_eV 		= 27.21139
au_to_angtrom	= 0.529177211


seed_name		= 'wf1'



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
	whitespace = '\t'

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
	print('wrote '+seed_name+'_hr.dat'+' file')


def write_r_file(seed_name, nAt, rhopp ):
	whitespace = '\t'
	rhopp = sorted(rhopp)

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
	print('wrote '+seed_name+'_r.dat'+' file')

#
#
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#		MEPinterp CFG FILE
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def write_mepInterp_input(file_path,valence_bands, ax, ay, az, a0, mp_grid, seed_name, kubo_tol=1e-3, hw=0.0,eFermi=0.0, Tkelvin=0.0,eta_smearing=0.0,  plot_bands='F',	debug_mode='F', do_gauge_trafo='T'):
	with open(file_path+'input.cfg','w') as outfile:
		outfile.write('# input file for TB model from New J Physics 12, 053032 (2010)'+'\n')
		outfile.write('# generated on '+datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")+'\n')
		outfile.write('\n')
		#
		#
		outfile.write('[jobs]\n')
		outfile.write('    '	+	'plot_bands='		+	str(plot_bands)			+	'\n')
		outfile.write('    '	+	'debug_mode='		+	str(debug_mode)			+	'\n')		
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
		#outfile.write('    '	+	'use_interp_kpt= '	+	use_interp_kpt		+	'\n')
		#outfile.write('    '	+	'do_gauge_trafo= '	+	do_gauge_trafo		+	'\n')
		outfile.write('\n')
		#
		#
		outfile.write('[MEP]\n')
		outfile.write('    '	+	'valence_bands= '	+	str(valence_bands)	+	'\n')
		#
		#
		outfile.write('[Kubo]\n')
		outfile.write('    '	+	'kuboTol= '			+	str(kubo_tol)		+	'\n')
		outfile.write('    '	+	'hw= '				+	str(hw)				+	'\n')
		outfile.write('    '	+	'eFermi= '			+	str(eFermi)			+	'\n')
		outfile.write('    '	+	'Tkelvin= '			+	str(Tkelvin)		+	'\n')
		outfile.write('    '	+	'eta_smearing= '	+	str(eta_smearing)	+	'\n')


		print('wrote '+file_path+'input.cfg')





#
#
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#		PUBLIC FUNCTION
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def write_souza_tb_input(root_dir, phi_para, valence_bands, mp_grid , kubo_tol=1e-3, hw=0.0, eFermi=0.0, Tkelvin=0.0, eta_smearing=0.0, plot_bands='F', debug_mode='F' ,do_gauge_trafo='T'):
	target_dir_name	= 'w90files'
	target_path		= root_dir+'/'+target_dir_name
	#
	#		get the tight binding basis
	nWfs, nrpts, tHopp, rHopp	=	get_souza_tb(phi_para)
	#
	#
	#write them to file
	os.mkdir(target_path)
	write_hr_file(	target_path+'/'+seed_name,	nWfs, nrpts, tHopp )
	write_r_file(	target_path+'/'+seed_name, 	nWfs,		 rHopp )

	# now write the input files for postw90 & mepInterp
	write_postw90_input(target_path, seed_name, valence_bands, mp_grid,  hw, eFermi, Tkelvin, eta_smearing	)

	#lattice setup
	ax 				= np.zeros(3)
	ay 				= np.zeros(3)
	az 				= np.zeros(3)
	ax[0]			= 2.0
	ay[1]			= 2.0
	az[2]			= 2.0
	a0				= 1.0
	#	write mepInterp
	write_mepInterp_input( root_dir+'/',valence_bands, ax, ay, az, a0, mp_grid, seed_name,kubo_tol, hw,eFermi, Tkelvin,eta_smearing, plot_bands, debug_mode, do_gauge_trafo)




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

	write_souza_tb_input(root_dir, phi_para, valence_bands, mp_grid)


#uncomment next line to generate some test input
test()


