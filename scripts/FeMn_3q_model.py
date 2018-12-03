import numpy as np
import 	datetime
#
#
#	parameter to convert from degrees to radiants
conv_deg_to_rad	=	np.pi / 180.
#
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#		CLASS HOLDING THE PARAMETERS																											   |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class FeMn_3q_model:
	#
	#
	def __init__(self, fpath,	verbose=False):
		self.nWfs		=	8
		self.nrpts		=	1
		self.verbose	=	verbose
		#
		self.thconv		=	[]
		self.phconv		=	[] 
		self.wf_pos		=	[]
		self.R_nn_lst	=	[]
		#
		self.v_print("	~~~~~~~~~~~~3D FCC - FeMn - 3q state~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		#
		#	consider only hopping with nearest neigbours
		#self.R_nn_lst.append(	[-1,  0,  0]	)
		#self.R_nn_lst.append(	[ 0, -1,  0]	)
		#self.R_nn_lst.append(	[ 0,  0, -1]	)
		self.R_nn_lst.append(	[ 0,  0,  0]	)
		#self.R_nn_lst.append(	[ 1,  0,  0]	)
		#self.R_nn_lst.append(	[ 0,  1,  0]	)
		#self.R_nn_lst.append(	[ 0,  0,  1]	)
		#
		#	real space positions
		rx_at	=	[.0,	.5,		.5,		.0]
		ry_at	=	[.0,	.5,		.0,		.5]
		rz_at	=	[.0,	.0,		.5,		.5]

		#	for spin up and down
		rx_wf	=	rx_at	+	rx_at
		ry_wf	=	ry_at	+	ry_at
		rz_wf	=	rz_at	+	rz_at

		print("[FeMn_3q_model/init]: rx_wf=",rx_wf)
		#
		for wf in range(self.nWfs):
			#position op. is real
			self.wf_pos.append(		[	rx_wf[wf],.0,		ry_wf[wf],.0,		rz_wf[wf],.0	]		)
		#
		#
		#	spin structure
		self.v_print(		'		  S4 ------ S2  '			)
		self.v_print(		'		   \        /   '			)
		self.v_print(		'		    \  S1  /    '			)
		self.v_print(		'		y    \    /     '			)
		self.v_print(		'		|     \  /      '			)
		self.v_print(		'		|      S3       '			)
		self.v_print(		'		O---x           '			)
		#
		with open(fpath) as f_in:
			for idx, line in enumerate(f_in):
				#	
				#	IDENTIFY INDIVIDUAL LINES 
				#
				if idx	== 0:
					self.nKx, self.nKy, self.nKz				=		np.fromstring(	line,		dtype=float, count=3, sep=" ")			
				#-------------------------------------------------------------------------------------
				#
				elif idx == 1:
					self.t1, self.t2, self.lmbda				=		np.fromstring(	line,		dtype=float, count=3, sep=" ")
				#-------------------------------------------------------------------------------------
				#
				elif idx >=2 and idx<=5:
					phi, theta									= 		np.fromstring(	line,		dtype=float, count=2, sep=" ")		
					self.phconv.append(		phi			* conv_deg_to_rad	)
					self.thconv.append(		theta		* conv_deg_to_rad	)
				#-------------------------------------------------------------------------------------
				#-------------------------------------------------------------------------------------
				#-------------------------------------------------------------------------------------
		self.v_print(	"	[FeMn_3q_model]:	initalized new param set from "+fpath			)
		self.v_print(	"	[FeMn_3q_model]:	 summary"										)
		self.v_print(	"		t1	    =	" + str(	self.t1		)	+ str("(eV)")		)
		self.v_print(	"		t2	    =	" + str(	self.t2		)	+ str("(eV)")		)
		self.v_print(	"		lmbda	=	" + str(	self.lmbda	)	+ str("(eV)")		)
		self.v_print(	" 		# spin		|	phi(rad) 		| theta(rad) "					)
		self.v_print(	"		------------------------------------------------------"			)
		for  spin in range(4):
			self.v_print("		  "+		str(spin+1)	+ "		|	"	+	str(self.phconv[spin])	+	"		|	"	+	str(self.thconv[spin])	)
		self.v_print(	"		------------------------------------------------------"			)
		self.v_print(	"		------------------------------------------------------"			)
	#
	#
	#	----
	def v_print(self,out_string):
		#print only if verbose option is true
		if self.verbose:	print(out_string)
	







	def tHopp_fill_zeros(self, tHopp):
		#	FILL REST WITH ZERO
		#
		# creat list with already added elements
		tHopp_exist = []
		exist_cnt	= 0
		for t in tHopp:
			tHopp_exist.append(	[	int(t[0]),int(t[1]),int(t[2]),int(t[3]),int(t[4])	]	 )
			exist_cnt	= exist_cnt + 1
		if exist_cnt != 80:
			print("WARNING the tHopp_exist list is wrong")
		#
		#only if not in already added list append value
		zero_cnt = 0
		match_cnt = 0
		t_matches = []
		for R_nn in self.R_nn_lst:
			for m in range(self.nWfs):
				for n in range(self.nWfs):
					if [ int(R_nn[0]), int(R_nn[1]), int(R_nn[2]), int(m+1), int(n+1)] in tHopp_exist:
						t_matches.append([ int(R_nn[0]), int(R_nn[1]), int(R_nn[2]), int(m+1), int(n+1)])
						match_cnt = match_cnt + 1				
					else:			
						zero_cnt	= zero_cnt + 1
						tHopp.append(		[ R_nn[0], R_nn[1], R_nn[2], m+1, n+1, 0.0, 0.0	]		)
		return tHopp















	#
	#
	#	----
	def setup_Ham(self):
		tHopp	=	[]
		#tHopp.append(	[ 0, 0, 0,			1,	2,			np.real(hopping[x][at1])	,   np.imag(hopping[x][at1])	])
		#
		#
		#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		#      EXCHANGE																|
		#---------------------------------------------------------------------------
		sintheta	= []
		costheta	= []
		phasphi		= []
		for i in range(4):
			sintheta.append(						np.sin( self.thconv[i] )							)
			costheta.append(						np.cos(	self.thconv[i] )							)
			phasphi.append(		np.complex128(	np.cos(self.phconv[i])	 - 1j * np.sin(self.phconv[i]))	) 	
			#
			z_upDw	=	self.lmbda 	* sintheta[i]	*	phasphi[i]
			re_upDw	=	np.real(	z_upDw	)
			im_upDw	=	np.imag(	z_upDw	)
			#       do i=1,4 ! exchange
			#        ham(i,i) = lambda*costheta(i)
			#        ham(i+4,i+4) = -lambda*costheta(i)
			#        ham(i,i+4) = lambda*phasphi(i)*sintheta(i)
			#       enddo
			tHopp.append(	[	0, 0, 0, 		i+1, i+1, 		self.lmbda 	* costheta[i]						, 	.0			]				)
			tHopp.append(	[ 	0, 0, 0,		i+5, i+5, 	-	self.lmbda	* costheta[i] 						, 	.0			]				)
			tHopp.append(	[	0, 0, 0,		i+1, i+5,			re_upDw										, 	im_upDw		]				)
			tHopp.append(	[	0, 0, 0, 		i+5, i+1,			re_upDw										, - im_upDw		]				)
		#
		#
		#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		#      HOPPING																|
		#---------------------------------------------------------------------------
		
		#! intra-layer hopping
       	#ht1(1) = -2*t*cos(pi*kx-3./2.*pi*ky)
       	#ht1(2) = -2*t*cos(2.*pi*kx)
       	#ht1(3) = -2*t*cos(pi*kx+3./2.*pi*ky)
       	

       	#! inter-layer hopping
       	#ht2(1) =-2*t2*cos(-pi*kx-pi/2.*ky-2.*pi*kz)
       	#ht2(2) =-2*t2*cos(          pi*ky-2.*pi*kz)
       	#ht2(3) =-2*t2*cos( pi*kx-pi/2.*ky-2.*pi*kz)

		#do i=0,4,4 ! intra- and inter-layer hopping
        #	ham(1+i,2+i)=ht1(1)+ht2(1)
        #	ham(1+i,3+i)=ht1(2)+ht2(2)
        #	ham(1+i,4+i)=ht1(3)+ht2(3)
        #	ham(2+i,3+i)=ht1(3)+ht2(3)
        #	ham(2+i,4+i)=ht1(2)+ht2(2)
        #	ham(3+i,4+i)=ht1(1)+ht2(1)
       	#enddo

		#
		#
		#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		#      FILL REST WITH ZEROS													|
		#---------------------------------------------------------------------------
		tHopp	=	self.tHopp_fill_zeros(tHopp)	
		return tHopp	
	










































	#
	#
	#	----
	def setup_Pos(self):
		rHopp		=	[]
		#
		self.v_print(	"		[FeMn_3q_model/setup_Pos]:	start position operator setup!"		)
		self.v_print(	"		[FeMn_3q_model/setup_Pos]: got the following nn cells: "		)
		self.v_print(	"		 n	|	R_nn")
		self.v_print(	"		-----------------------------------")
		for idx,R_nn in enumerate(self.R_nn_lst):
			self.v_print(" 		"+str(idx)+"  |  "+str(R_nn))
		self.v_print(	"		-----------------------------------")

		#
		for m in range(self.nWfs):
			for n in range(self.nWfs):
				for idx,R_nn in enumerate(self.R_nn_lst):
					at_pos 	=	[	.0,.0,		.0,.0,		.0,.0]
					#	R_nn is diagonal
					if R_nn == [0,0,0]:		
						if m == n:
							at_pos	=	self.wf_pos[m]		
					rHopp.append(	[	R_nn[0],R_nn[1],R_nn[2], 	m+1, n+1, 		at_pos[0], at_pos[1], at_pos[2], at_pos[3], at_pos[4], at_pos[5]		])
		return rHopp
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






#
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#		INTERFACE																																   |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_FeMn3q_tb(fpath, verbose=False):
	print("[get_FeMn3q_tb]:	started..  FeMn3q model setup at " +			datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"))
	#
	#	READ INPUT FILE
	sys_params		=	FeMn_3q_model(fpath, verbose)
	#
	#
	#	HOPPING
	tHopp			=	sys_params.setup_Ham()
	#
	#	POSITION
	rHopp			=	sys_params.setup_Pos()
	#
	#	LATTICE
	a0 	= 1.0
	ax	= np.zeros(3)
	ay 	= np.zeros(3)
	az	= np.zeros(3)
	ax[0]	=	1.0
	ay[1]	=	1.0
	az[2]	=	1.0	
	#
	#
	print("[get_FeMn3q_tb]:	..finished FeMn3q model setup (TODO: PROPER LATTICE SETUP !) at " +	datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"))
	#
	return sys_params.nWfs, sys_params.nrpts, tHopp, rHopp, ax, ay, az, a0


#
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#		TESTING																																	   |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def test(verbose=True):
	fpath 	= '/Users/merte/bin/tbmodel_3qstate/inp_params_3q'
	#
	#
	print("*")
	print("*")
	print("*")
	print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|")
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TEST_OF_FeMn_3q_CLASS~~~~~~~~~~|")
	print("")
	print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|")
	print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|")
	if verbose:	print("[FeMn_3q_model/test]:	verbose mode active")
	print("")
	print("")
	print("")
	#
	#
	nWfs, nrpts, tHopp, rHopp, ax, ay, az, a0	=	get_FeMn3q_tb(fpath, verbose)
	#
	#
	#
	print("")
	print("")
	print("")
	print("*")
	print("*")
	print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|")
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SUMMARY~~~~~~~~~~~~~~~~~~~~~~~~|")
	print("")
	print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|")
	print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|")
	print("[FeMn_3q_model/test]:			nWfs ="	+	str(	nWfs			))
	print("[FeMn_3q_model/test]:			nrpts="	+	str(	nrpts			))
	print("	"	)
	print("LATTICE")
	print("ax=	",a0*ax,"	(arb. units)")
	print("ay=	",a0*ay,"	(arb. units)")
	print("az=	",a0*az,"	(arb. units)")


	print("[FeMn_3q_model/test]:			tHopp (len="+str(len(tHopp))+"): [eV]")
	print("")
	print("	R 	|	n m  |	t^R_nm (eV)")
	print("----------------------------------------------")
	for t_mn in tHopp:
		print(t_mn[0:3],"	| ",t_mn[3:5],"	|	",t_mn[-1])
	
	print("")
	print("[FeMn_3q_model/test]:			rHopp (len="+str(len(rHopp))+")")
	#for r_nm in rHopp:
	#	print(r_nm)
	print("~~~~~~~~~~~~~~~~~~~~~~~~")
	

#
#
#UNCOMMENT TO TEST THIS SCRIPT
#test(verbose	=	False)







