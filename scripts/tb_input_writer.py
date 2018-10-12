import numpy as np
import datetime
import os
from	shutil import rmtree
from postw90_in_writer import postw90_job


au_to_eV 		= 27.21139
au_to_angtrom	= 0.529177211


seed_name		= 'wf1'

#--------------------MEPinterp input paras--------------------------------------------------------------------------------
#write_mepInterp_input(plot_bands, a1, a2, a3, a0, mp_grid, seed_name, use_interp_kpt, do_gauge_trafo, valence_bands)




#
#----------------------FUNCTIONS-------------------------------------------------------------------------------------------

def convert_phase_to_complex(angles):
	#computes e^{i*phi}
	tmp = []
	for angle in angles:
		tmp.append(	np.exp(	1j * np.pi * angle)		)
		#tmp.append(		2.0*np.cos(np.pi*angle)			)
	return tmp




def get_nn_dist(atom_pos):
	#search for smallest distance between different elements within atom_pos
	nn_dist = 0.0
	#
	for at, at_pos in enumerate(atom_pos):
		for nn, nn_pos in enumerate(atom_pos):
			if at is not nn:
				#
				tmp = np.linalg.norm( nn_pos - at_pos)
				if nn==2 and at==1:
					nn_dist = tmp
				elif tmp < nn_dist:
					nn_dist = tmp
	return nn_dist




def int_to_dim_string(x):
	#give string denoting the coordinate axis
	if x== 0:
		dim_string = 'x'
	elif x == 1:
		dim_string = 'y'
	elif x == 2:
		dim_string = 'z'
	else:
		dim_string = 'this shit is of higher dimension'

	return dim_string



def get_at_pos(a_latt, rel_atom_pos):
	rel_vec = np.array(rel_atom_pos)
	pos =  a_latt.dot(rel_vec)
	#
	cmplx_pos = np.zeros(6)
	cmplx_pos[0]	= pos[0]
	cmplx_pos[2]	= pos[1]
	cmplx_pos[4]	= pos[2]

	return cmplx_pos


def get_hopp_list(nAt, rel_atom_pos, onsite, phi):
	#setup the hopping list
	thopp	= []

	#add the onsite energies to the hopping list
	for at, en in enumerate(onsite):
		thopp.append([0,0,0,at+1,at+1,en,.0])


	#add nearest neigbour hopping to hopping list
	nn_dist = get_nn_dist(rel_atom_pos)
	#print('determined (relative) nn_dist:'+str(nn_dist))
	#
	nn_list = [[None for x in range(nAt)] for y in range(nAt)]
	for at in range(nAt):
		for nn, nn_pos in enumerate(rel_atom_pos):
			b_nn	= nn_pos - rel_atom_pos[at]
			nn_list[at][nn]	= '0'
			#found nearest neighbour
			if abs(np.linalg.norm(b_nn)-nn_dist) < 1e-8:
				#determine R_nn
				R_nn = [0,0,0]
				for x in range(3):
					# nn is in home unit cell (R_nn=0)
					if b_nn[x] > 1e-8:
						b_dim 	= x
						nn_list[at][nn]	= int_to_dim_string(b_dim)
					# nn is in unit cell to the right
					elif b_nn[x] < - 1e-8:
						b_dim 	= x
						R_nn[x]	= int(	np.sign(b_nn[x])		)
						nn_list[at][nn]	= int_to_dim_string(b_dim)


				#print('at='+str(at+1)+' nn='+str(nn+1)+' b_nn('+nn_list[at][nn]+')='+str(b_nn)+' R_nn='+str(R_nn)+' phi='+str(phi[b_dim][at]))
				thopp.append(	[		R_nn[0], R_nn[1], R_nn[2], at+1, nn+1, np.real(phi[b_dim][at]), +np.imag(phi[b_dim][at])		]	)
				thopp.append(	[		R_nn[0], R_nn[1], R_nn[2], nn+1, at+1, np.real(phi[b_dim][at]), -np.imag(phi[b_dim][at])		]	)

				

	#for at,nn in enumerate(nn_list):
	#	print('at='+str(at)+' nn='+str(nn)	)


	#get a list of all R_hopp vectors used
	R_hopp = []
	for hopp in thopp:
		if [hopp[0],hopp[1],hopp[2]] in R_hopp:
			x = x
		else:
			R_hopp.append([hopp[0],hopp[1],hopp[2]])

	nrpts= len(R_hopp)

	#list with already added elements
	thopp_helper = []
	for t in thopp:
		thopp_helper.append(	[	t[0],t[1],t[2],t[3],t[4]	]	 )

	#only if not in already added list append value
	for R_nn in R_hopp:
		for m in range(nAt):
			for n in range(nAt):
				if [ R_nn[0], R_nn[1], R_nn[2], m+1, n+1] not in thopp_helper:
					thopp.append(		[ R_nn[0], R_nn[1], R_nn[2], m+1, n+1, 0.0, 0.0	]		)

	#sort targets
	return nrpts, thopp, R_hopp



def get_pos_list(nAt, rel_atom_pos, a_latt, R_nn_list):
	rhopp = []

	for R_nn in R_nn_list:
		for m in range(nAt):
			for n in range(nAt):
				if abs(np.linalg.norm(R_nn)) < 1e-9 and m==n:
					at_pos = get_at_pos(a_latt, rel_atom_pos[m])
				else:
					at_pos = np.zeros(6)
				rhopp.append(	[		R_nn[0],R_nn[1],R_nn[2], m+1, n+1, at_pos[0], at_pos[1], at_pos[2], at_pos[3], at_pos[4], at_pos[5]			]		)
	return rhopp



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




def write_mepInterp_input(file_path,valence_bands, ax, ay, az, a0, mp_grid, seed_name, kubo_tol=1e-3, hw=0.0,eFermi=0.0, Tkelvin=0.0,eta_smearing=0.0,  plot_bands='F',	debug_mode='F', do_gauge_trafo='T'):
	with open(file_path+'input.txt','w') as outfile:
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


		print('wrote '+file_path+'input.txt')



def write_postw90_input(w90_dir, seed_name, val_bands, mp_grid, hw=0.0, eFermi=0.0, Tkelvin=0.0, eta_smearing=0.0):
	pw90_job 	=	postw90_job(w90_dir, seed_name, val_bands, mp_grid,  hw, eFermi, Tkelvin, eta_smearing)
	pw90_job.write_win_file()
	#pw90_job.run()
#
#----------------------BODY------------------------------------------------------------------------------------------------






def write_souza_tb_input(root_dir, phi_para, valence_bands, mp_grid , kubo_tol=1e-3, hw=0.0, eFermi=0.0, Tkelvin=0.0, eta_smearing=0.0, plot_bands='F', debug_mode='F' ,do_gauge_trafo='T'):
	target_dir_name	= 'w90files'
	target_path		= root_dir+'/'+target_dir_name


	#----------------------TB LATTICE---------------------------------------------------------------------------------------
	ax 				= np.zeros(3)
	ay 				= np.zeros(3)
	az 				= np.zeros(3)
	ax[0]			= 2.0
	ay[1]			= 2.0
	az[2]			= 2.0
	a0				= 1.0
	a_latt 			= np.array(
					[
						[ax[0], ax[1], ax[2]],
						[ay[0],	ay[1], ay[2]],
						[az[0], az[1], az[2]],
					])
	a_latt			= a_latt * au_to_angtrom

	#----------------------TB ATOMS---------------------------------------------------------------------------------------
	nAt				= 8
	rel_atom_pos 	= np.array(
					[	[.0,.0,.0],
						[.5,.0,.0],
						[.5,.5,.0],
						[.0,.5,.0],
						[.0,.0,.5],
						[.5,.0,.5],
						[.5,.5,.5],
						[.0,.5,.5]
					])

	if valence_bands > nAt:
		print('WARNING: more valence_bands then total number of bands('+str(nAt)+') will set valence_bands to '+str(nAt))
		valence_bands = nAt

	#----------------------TB PARAMETERS---------------------------------------------------------------------------------------
	onsite	=	np.array(	[-6.5, 0.9, 1.4, 1.2, -6.0, 1.5, 0.8, 1.2]			)
	phi_x	=	np.array(	[phi_para, 	1.3, 0.8, 0.3, 1.4, 0.6, 0.8, 1.9]		)
	phi_y	=	np.array(	[	0.5,	0.2, 1.4, 1.9, 0.8, 1.7, 0.6, 0.3]		)
	phi_z	=	np.array(	[	1.7,	0.5, 0.6, 1.0, 0.3, 0.7, 1.2, 1.4]		)




	#prepare hopping values
	phi_x	= convert_phase_to_complex(phi_x)
	phi_y	= convert_phase_to_complex(phi_y)
	phi_z	= convert_phase_to_complex(phi_z)
	


	#
	phi		= []
	phi.append(phi_x)
	phi.append(phi_y)
	phi.append(phi_z)

	onsite	= onsite	* au_to_eV
	phi 	= np.array(phi)	      * au_to_eV
	#print('hoppings (phi='+str(phi_para)+'):	')
	#print(phi)
	#get the lists
	nrpts, thopp, R_nn_list = get_hopp_list(nAt, rel_atom_pos, onsite, phi)
	rhopp	= get_pos_list(nAt, rel_atom_pos, a_latt, R_nn_list)




	#write them to file
	os.mkdir(target_path)
	write_hr_file(	target_path+'/'+seed_name,	nAt, nrpts, thopp )
	write_r_file(	target_path+'/'+seed_name, 	nAt,		rhopp )

	# now write the input files for postw90 & mepInterp
	write_postw90_input(target_path, seed_name, valence_bands, mp_grid,  hw, eFermi, Tkelvin, eta_smearing	)


	write_mepInterp_input( root_dir+'/',valence_bands, ax, ay, az, a0, mp_grid, seed_name,kubo_tol, hw,eFermi, Tkelvin,eta_smearing, plot_bands, debug_mode, do_gauge_trafo)


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
#test()


