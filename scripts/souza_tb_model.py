import numpy as np

#---- PARAMETERS -----------------------
#
au_to_eV 		= 27.21139
au_to_angtrom	= 0.529177211
#
x	= 0
y 	= 1
z 	= 2
#
at1 = 0
at2	= 1
at3 = 2
at4 = 3
at5 = 4
at6 = 5
at7 = 6
at8 = 7
#******************************************




def convert_phase_to_complex(angles):
	#computes e^{i*phi}
	hopp = []
	for angle in angles:
		hopp.append(	np.exp(	1j * np.pi * angle)		)
		#tmp.append(		2.0*np.cos(np.pi*angle)			)
	return hopp


def prepare_hoppings(phi_x, phi_y, phi_z):
	#	get 		exp(i	Phi)
	phi_x	= convert_phase_to_complex(phi_x)
	phi_y	= convert_phase_to_complex(phi_y)
	phi_z	= convert_phase_to_complex(phi_z)
	#
	#	fill target array
	phi		= []
	phi.append(phi_x)
	phi.append(phi_y)
	phi.append(phi_z)
	#
	return np.array(phi)

def set_units(onsite, phi):
	onsite	= 	onsite			* au_to_eV
	phi 	= 	phi			    * au_to_eV
	#
	return onsite, phi



#	HOPPING HELPERS:


	#	*********************************************************
	#		schematic unit cell 
	#		#5-#8 indicate a second layer in the unit cell along z 
	#	-----
	#
	#
	#									#8 - #7
	#								  /	|	  |
	#								/	#5 - #6
	#							  /	    /    /
	#							/	  /	    /  
	#		    			 /	    /     /
	#						#4 - #3     /
	#	 y					|	  |   /
	#	 ^	 				#1 - #2 /
	#	 |  
	#	 |  
	#	 |---> x
	#	------------------------------------------------------------



def get_x_hopp(hopping):
	tHopp	=	[]
	#	
	#			X-HOPPINGS
	#
	#	atom 1 - 2 		Rx
	tHopp.append(	[-1, 0, 0,			1,	2,			np.real(hopping[x][at1])	,	np.imag(hopping[x][at1])	])
	tHopp.append(	[ 0, 0, 0,			2,	1,			np.real(hopping[x][at2])	,	np.imag(hopping[x][at2])	])
	tHopp.append(	[+1, 0, 0,			1,	2,			np.real(hopping[x][at1])	,	np.imag(hopping[x][at1])	])
	#	c.c.:
	tHopp.append(	[-1, 0, 0,			2,	1,			np.real(hopping[x][at1])	, - np.imag(hopping[x][at1])	])
	tHopp.append(	[ 0, 0, 0,			1,	2,			np.real(hopping[x][at2])	, - np.imag(hopping[x][at2])	])
	tHopp.append(	[+1, 0, 0,			2,	1,			np.real(hopping[x][at1])	, - np.imag(hopping[x][at1])	])
	#
	#	atom 4 - 3		Rx
	tHopp.append(	[-1, 0, 0,			4,	3,			np.real(hopping[x][at4])	,	np.imag(hopping[x][at4])	])
	tHopp.append(	[ 0, 0, 0,			3,	4,			np.real(hopping[x][at3])	,	np.imag(hopping[x][at3])	])
	tHopp.append(	[+1, 0, 0,			4,	3,			np.real(hopping[x][at4])	,	np.imag(hopping[x][at4])	])
	#	c.c.:
	tHopp.append(	[-1, 0, 0,			3,	4,			np.real(hopping[x][at4])	, - np.imag(hopping[x][at4])	])
	tHopp.append(	[ 0, 0, 0,			4,	3,			np.real(hopping[x][at3])	, - np.imag(hopping[x][at3])	])
	tHopp.append(	[+1, 0, 0,			3,	4,			np.real(hopping[x][at4])	, - np.imag(hopping[x][at4])	])
	#
	#
	#	atom 5 - 6 		Rx
	tHopp.append(	[-1, 0, 0,			5,	6,			np.real(hopping[x][at5])	,	np.imag(hopping[x][at5])	])
	tHopp.append(	[ 0, 0, 0,			6,	5,			np.real(hopping[x][at6])	,	np.imag(hopping[x][at6])	])
	tHopp.append(	[+1, 0, 0,			5,	6,			np.real(hopping[x][at5])	,	np.imag(hopping[x][at5])	])
	#	c.c.:
	tHopp.append(	[-1, 0, 0,			6,	5,			np.real(hopping[x][at5])	, - np.imag(hopping[x][at5])	])
	tHopp.append(	[ 0, 0, 0,			5,	6,			np.real(hopping[x][at6])	, - np.imag(hopping[x][at6])	])
	tHopp.append(	[+1, 0, 0,			6,	5,			np.real(hopping[x][at5])	, - np.imag(hopping[x][at5])	])
	#
	#
	#	atom 8 - 7		Rx
	tHopp.append(	[-1, 0, 0,			8,	7,			np.real(hopping[x][at8])	,	np.imag(hopping[x][at8])	])
	tHopp.append(	[ 0, 0, 0,			7,	8,			np.real(hopping[x][at7])	,	np.imag(hopping[x][at7])	])
	tHopp.append(	[+1, 0, 0,			8,	7,			np.real(hopping[x][at8])	,	np.imag(hopping[x][at8])	])
	#	c.c.:
	tHopp.append(	[-1, 0, 0,			7,	8,			np.real(hopping[x][at8])	, - np.imag(hopping[x][at8])	])
	tHopp.append(	[ 0, 0, 0,			8,	7,			np.real(hopping[x][at7])	, - np.imag(hopping[x][at7])	])
	tHopp.append(	[+1, 0, 0,			7,	8,			np.real(hopping[x][at8])	, - np.imag(hopping[x][at8])	])
	#
	return tHopp 


def get_y_hopp(hopping):
	tHopp	=	[]	
	#			Y-HOPPINGS
	#
	#	atom 1 - 4 		Ry
	tHopp.append(	[0, -1, 0,			1,	4,			np.real(hopping[y][at1])	,	np.imag(hopping[y][at1])	])
	tHopp.append(	[0,  0, 0,			4,	1,			np.real(hopping[y][at4])	,	np.imag(hopping[y][at4])	])
	tHopp.append(	[0, +1, 0,			1,	4,			np.real(hopping[y][at1])	,	np.imag(hopping[y][at1])	])
	#	c.c.:
	tHopp.append(	[0, -1, 0,			4,	1,			np.real(hopping[y][at1])	, - np.imag(hopping[y][at1])	])
	tHopp.append(	[0,  0, 0,			1,	4,			np.real(hopping[y][at4])	, - np.imag(hopping[y][at4])	])
	tHopp.append(	[0, +1, 0,			4,	1,			np.real(hopping[y][at1])	, - np.imag(hopping[y][at1])	])
	#
	#
	#	atom 2 - 3 		Ry
	tHopp.append(	[0, -1, 0,			2,	3,			np.real(hopping[y][at2])	,	np.imag(hopping[y][at2])	])
	tHopp.append(	[0,  0, 0,			3,	2,			np.real(hopping[y][at3])	,	np.imag(hopping[y][at3])	])
	tHopp.append(	[0, +1, 0,			2,	3,			np.real(hopping[y][at2])	,	np.imag(hopping[y][at2])	])
	#	c.c.:
	tHopp.append(	[0, -1, 0,			3,	2,			np.real(hopping[y][at2])	, - np.imag(hopping[y][at2])	])
	tHopp.append(	[0,  0, 0,			2,	3,			np.real(hopping[y][at3])	, - np.imag(hopping[y][at3])	])
	tHopp.append(	[0, +1, 0,			3,	2,			np.real(hopping[y][at2])	, - np.imag(hopping[y][at2])	])
	#
	#
	#	atom 5 - 8 		Ry
	tHopp.append(	[0, -1, 0,			5,	8,			np.real(hopping[y][at5])	,	np.imag(hopping[y][at5])	])
	tHopp.append(	[0,  0, 0,			8,	5,			np.real(hopping[y][at8])	,	np.imag(hopping[y][at8])	])
	tHopp.append(	[0, +1, 0,			5,	4,			np.real(hopping[y][at5])	,	np.imag(hopping[y][at5])	])
	#	c.c.:
	tHopp.append(	[0, -1, 0,			8,	5,			np.real(hopping[y][at5])	, - np.imag(hopping[y][at5])	])
	tHopp.append(	[0,  0, 0,			5,	8,			np.real(hopping[y][at8])	, - np.imag(hopping[y][at8])	])
	tHopp.append(	[0, +1, 0,			8,	5,			np.real(hopping[y][at5])	, - np.imag(hopping[y][at5])	])
	#
	#
	#	atom 6 - 7 		Ry
	tHopp.append(	[0, -1, 0,			6,	7,			np.real(hopping[y][at6])	,	np.imag(hopping[y][at6])	])
	tHopp.append(	[0,  0, 0,			7,	6,			np.real(hopping[y][at7])	,	np.imag(hopping[y][at7])	])
	tHopp.append(	[0, +1, 0,			6,	7,			np.real(hopping[y][at6])	,	np.imag(hopping[y][at6])	])
	#	c.c.:
	tHopp.append(	[0, -1, 0,			7,	6,			np.real(hopping[y][at6])	, - np.imag(hopping[y][at6])	])
	tHopp.append(	[0,  0, 0,			6,	7,			np.real(hopping[y][at7])	, - np.imag(hopping[y][at7])	])
	tHopp.append(	[0, +1, 0,			7,	6,			np.real(hopping[y][at6])	, - np.imag(hopping[y][at6])	])
	#
	return tHopp



def get_z_hopp(hopping):
	tHopp 	=	[]
	#	
	#			Z-HOPPINGS
	#
	#	atom 1 - 5 		Rz
	tHopp.append(	[0, 0, -1,			1,	5,			np.real(hopping[z][at1])	,	np.imag(hopping[z][at1])	])
	tHopp.append(	[0, 0,  0,			5,	1,			np.real(hopping[z][at5])	,	np.imag(hopping[z][at5])	])
	tHopp.append(	[0, 0, +1,			1,	5,			np.real(hopping[z][at1])	,	np.imag(hopping[z][at1])	])
	#	c.c.:
	tHopp.append(	[0, 0, -1,			5,	1,			np.real(hopping[z][at1])	, - np.imag(hopping[z][at1])	])
	tHopp.append(	[0, 0,  0,			1,	5,			np.real(hopping[z][at5])	, - np.imag(hopping[z][at5])	])
	tHopp.append(	[0, 0, +1,			5,	1,			np.real(hopping[z][at1])	, - np.imag(hopping[z][at1])	])
	#
	#
	#	atom 2 - 6 		Rz
	tHopp.append(	[0, 0, -1,			2,	6,			np.real(hopping[z][at2])	,	np.imag(hopping[z][at2])	])
	tHopp.append(	[0, 0,  0,			6,	2,			np.real(hopping[z][at6])	,	np.imag(hopping[z][at6])	])
	tHopp.append(	[0, 0, +1,			2,	6,			np.real(hopping[z][at2])	,	np.imag(hopping[z][at2])	])
	#	c.c.:
	tHopp.append(	[0, 0, -1,			6,	2,			np.real(hopping[z][at2])	, - np.imag(hopping[z][at2])	])
	tHopp.append(	[0, 0,  0,			2,	6,			np.real(hopping[z][at6])	, - np.imag(hopping[z][at6])	])
	tHopp.append(	[0, 0, +1,			6,	2,			np.real(hopping[z][at2])	, - np.imag(hopping[z][at2])	])
	#
	#
	#	atom 3 - 7 		Rz
	tHopp.append(	[0, 0, -1,			3,	7,			np.real(hopping[z][at3])	,	np.imag(hopping[z][at3])	])
	tHopp.append(	[0, 0,  0,			7,	3,			np.real(hopping[z][at7])	,	np.imag(hopping[z][at7])	])
	tHopp.append(	[0, 0, +1,			3,	7,			np.real(hopping[z][at3])	,	np.imag(hopping[z][at3])	])
	#	c.c.:
	tHopp.append(	[0, 0, -1,			7,	3,			np.real(hopping[z][at3])	, - np.imag(hopping[z][at3])	])
	tHopp.append(	[0, 0,  0,			3,	7,			np.real(hopping[z][at7])	, - np.imag(hopping[z][at7])	])
	tHopp.append(	[0, 0, +1,			7,	3,			np.real(hopping[z][at3])	, - np.imag(hopping[z][at3])	])
	#
	#
	#	atom 4 - 8 		Rz
	tHopp.append(	[0, 0, -1,			4,	8,			np.real(hopping[z][at4])	,	np.imag(hopping[z][at4])	])
	tHopp.append(	[0, 0,  0,			8,	4,			np.real(hopping[z][at8])	,	np.imag(hopping[z][at8])	])
	tHopp.append(	[0, 0, +1,			4,	8,			np.real(hopping[z][at4])	,	np.imag(hopping[z][at4])	])
	#	c.c.:
	tHopp.append(	[0, 0, -1,			8,	4,			np.real(hopping[z][at4])	, - np.imag(hopping[z][at4])	])
	tHopp.append(	[0, 0,  0,			4,	8,			np.real(hopping[z][at8])	, - np.imag(hopping[z][at8])	])
	tHopp.append(	[0, 0, +1,			8,	4,			np.real(hopping[z][at4])	, - np.imag(hopping[z][at4])	])
	#
	return tHopp












def get_souza_tb(phi_para):
	#
	onsite			=	np.array(	[-6.5, 0.9, 1.4, 1.2, -6.0, 1.5, 0.8, 1.2]			)
	phi_x			=	np.array(	[phi_para, 	1.3, 0.8, 0.3, 1.4, 0.6, 0.8, 1.9]		)
	phi_y			=	np.array(	[	0.5,	0.2, 1.4, 1.9, 0.8, 1.7, 0.6, 0.3]		)
	phi_z			=	np.array(	[	1.7,	0.5, 0.6, 1.0, 0.3, 0.7, 1.2, 1.4]		)
	#
	#
	hopping 		=	prepare_hoppings(phi_x, phi_y, phi_z)
	onsite, hopping	=	set_units(onsite, hopping)
	#
	#	INPUT DEBUG
	if onsite.size is not 8:
		print("[get_souza_tb]:	onsite energy arrray has wrong size")
		stop 
	if hopping.size is not 24:
		print("[get_souza_tb]:	hopping energy array has wrong size")
		stop
	#
	print('')
	print('')
	print('')
	print(' ~~~~~~~~~~~~~~ CUBIC TIGHT BINDING MODEL - NON INVERSION SYMMETRIC ~~~~~~~~~~~~~~~~~~~~')
	print('*')
	print('*')
	print('*')
	#
	#	SET MODEL PARAS
	nWfs	=	8
	nrpts	=	7
	#
	#
	#	INIT CONTAINERS
	R_hopp	=	[]
	#	0
	R_hopp.append(	[ 0,0,0]		)	# Home unit cell
	#	x
	R_hopp.append(	[-1,0,0]		)	# - X
	R_hopp.append(	[+1,0,0]		)	# + X
	#	x
	R_hopp.append(	[0,-1,0]		)	# - Y
	R_hopp.append(	[0,+1,0]		)	# + Y
	#	z
	R_hopp.append(	[0,0,-1]		)	# - Z
	R_hopp.append(	[0,0,+1]		)	# + Z
	#

	#
	#	ONSITE ENERGIES
	en_0	= []
	for at, en in enumerate(onsite):
		en_0.append(		[	0, 0, 0,			at+1, at+1,			en,	.0]			)
	
	#	HOPPING
	tHopp_X		=	get_x_hopp(hopping)
	tHopp_Y 	=	get_y_hopp(hopping)
	tHopp_Z		=	get_z_hopp(hopping)
	#
	#	COLLECT
	tHopp	=	np.concatenate(		(en_0,		tHopp_X),		axis=0)
	tHopp	=	np.concatenate(		(tHopp,		tHopp_Y),		axis=0)
	tHopp	=	np.concatenate(		(tHopp,		tHopp_Z),		axis=0)
	#
	#
	#	TEST IF ALL HOPPINGS WERE COLLECTED 
	exp_non_zero_hopp	=	8 + 3 * 24
	if int(tHopp.size/7) is not exp_non_zero_hopp: 			# 8 onsite energies, 24 hoppings per spatial dimension
		print("tHopp contains: "+str(tHopp.size/7)+" non zero values, expected 80")
	else:
		print("tHopp has "+str(exp_non_zero_hopp) +	" non zero entries")
	#
	#
	#	FILL REST WITH ZERO
	#
	# creat list with already added elements
	tHopp_exist = []
	exist_cnt	= 0
	for t in tHopp:
		tHopp_exist.append(	[	int(t[0]),int(t[1]),int(t[2]),int(t[3]),int(t[4])	]	 )
		exist_cnt	= exist_cnt + 1
	if exist_cnt is not 80:
		print("the tHopp_exist list is wrong")
	#
	#only if not in already added list append value
	zero_cnt = 0
	match_cnt = 0
	t_matches = []
	for R_nn in R_hopp:
		for m in range(nWfs):
			for n in range(nWfs):
				if [ int(R_nn[0]), int(R_nn[1]), int(R_nn[2]), int(m+1), int(n+1)] in tHopp_exist:
					t_matches.append([ int(R_nn[0]), int(R_nn[1]), int(R_nn[2]), int(m+1), int(n+1)])
					match_cnt = match_cnt + 1				
				else:			
					zero_cnt	= zero_cnt + 1

				#	zero_hopp	= np.array(	 [	[  R_nn[0], R_nn[1], R_nn[2], m+1, n+1, 0.0, 0.0 ]	 ]	)
				#	tHopp		= np.concatenate(	(tHopp, zero_hopp),		axis=0)	
					
									#tHopp = np.append( tHopp, np.array([  [R_nn[0], R_nn[1], R_nn[2], m+1, n+1, 0.0, 0.0]			])	)
	#				#tHopp.append(		[ R_nn[0], R_nn[1], R_nn[2], m+1, n+1, 0.0, 0.0	]		)
	

	if zero_cnt is not  nrpts*nWfs**2-exp_non_zero_hopp:
		print("detected wrong zero count, got "+str(zero_cnt)+" expected "+str(nrpts*nWfs**2-exp_non_zero_hopp))
		print("counted "+str(match_cnt)+" matches expected 80")
		print("t_matches:")
		for match in t_matches:
			print(str(match)+"	interpreted as non zero")
	else:
		print('added '+str(zero_cnt)+" zeros to the hopping matrix")

	if tHopp.size is not 7*nWfs**2 * nrpts:
		print("WARNING tHopp.size="+str(tHopp.size/7)+" expected "+str(nWfs**2 * nrpts))

	#
	#	POSITION



	return nWfs, nrpts, tHopp, R_hopp


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#
#
#
#

#
#
def test():
	nWfs, nrpts, tHopp, R_hopp 	=	get_souza_tb(0.0)
	#
	#print('tHopp:')
	#for t in tHopp:
	#	print(" "+str(int(t[0]))+ " "+str(int(t[1]))+ " "+str(int(t[2]))+" 		"+str(int(t[3]))+" "+str(int(t[4]))+" 		"+str(t[5])+	" "+str(t[6])	)
	#
	print('nWfs= '+str(nWfs))
	print('nrpts= '+str(nrpts))
	print('R_hopp= '+str(R_hopp))




test()
