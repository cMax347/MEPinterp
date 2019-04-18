import numpy as np
from scipy.optimize import minimize_scalar


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#	original MATLAB code from Wanxiang:
#
#
#		 1>	syms x_sym latt_sym vol_sym   % use symbolic computation
#		 2>	latt_sym = [ 2*x_sym          0      z;
#		 3>	              -x_sym   sqrt(3)*x_sym z;
#		 4>	              -x_sym  -sqrt(3)*x_sym z];
#		 5>	vol_sym = dot(latt_sym(1,:),cross(latt_sym(2,:),latt_sym(3,:)));
#		 6>	eqn = vol_sym == vol;
#		 7>	x_sym = double(solve(eqn,x_sym));
#		 8>	x_sym(real(x_sym)<0|imag(x_sym)~=0)=[];
#		 9>	
#		10>	x = double(x_sym);
#		11>	latt1 = [ 2*x        0    z;
#		12>	           -x   sqrt(3)*x z;
#		13>	           -x  -sqrt(3)*x z];
#		14>	latt = latt1;   % new lattice
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_vol(latt):
	vol	=	np.dot(latt[0,:],np.cross(latt[1,:],latt[2,:]))
	return vol	




class FeMn_latt_OPT:
	def print(self,msg):
		if self.verbose:
			print(msg)

	def __init__(self,strain, verbose=False):
		#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		#		ORIGINAL LATTICE																		|
		#																								|
		self.ndim 	= 3                                						   	# 	3D system			|
		self.a0		=	6.8597													#	lattice constant	|
		self.latt0 	= np.array([ 												#						|						
								[ 5.60092859,  0.00000000,  3.96045459],   	# 	original lattice		|
		         				[-2.80046430,  4.85054644,  3.96045459],   	# 	in unit of Bohr			|
			         			[-2.80046430, -4.85054644,  3.96045459]		#							|
			         			])												#						|
		self.z0 		= 	self.latt0[0,2]  	#  % original distance between adjacent (111) planes	|
		self.vol0		=	get_vol(self.latt0)								#							|
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|
		#	strained lattice
		self.strain		=	strain
		self.z1 		=	self.z0	* self.strain	
		#
		self.verbose	=	verbose
		#
		#	info
		self.print("")
		self.print("")
		self.print("^^^STRAINED LATTIC CALCULATOR^^^")
		self.print("\t target volume ="+str(self.vol0)	+'	a_0**3')
	#~~~~~~~~~~~~~~~~~
	#
	#
	def get_latt1(self,x):
		latt1	=	np.array(	[	[2.*x	,	0.				,		self.z1],
									[ - x	,	x*np.sqrt(3)	,		self.z1],
									[ - x	, -	x*np.sqrt(3)	,		self.z1]
								])
		return latt1
	#~~~~~~~~~~~~~~~~~
	#
	#
	def const_VOL(self,x):		
		#% solve the strained lattice, by keeping constant volume
		latt_x	=	self.get_latt1(x)
		vol1	=	get_vol(latt_x)
		return np.abs(vol1-self.vol0)
	#~~~~~~~~~~~~~~~~~
	#
	#
	def optimize_lattice(self,tol=1e-5):
		#print('[optimize_lattice]: proudly present the candidates:')
		self.print('\t [optimize_lattice]:  will use constant volume approximation')
		opt = minimize_scalar(		self.const_VOL, 	
									bounds=(0,self.a0),
									method='bounded', 
									tol=None,             
									options={'xatol': tol,'maxiter':5000,'disp': 0}
							)
		#
		if not opt.success:
			self.print('[optimize_lattice]: ERROR optimization failed: '+opt.message)
		opt_latt	=	self.get_latt1(		opt.x		)
		opt_vol		=	get_vol(		opt_latt		)
		#
		dV			=	np.abs(  self.vol0 - opt_vol )	
		if ( dV > tol*10.):
			self.print("[optimize_lattice]: WARNING  optimization did not converge dV="+str(dV))
		self.print('\t [optimize_lattice]:  optimized lattice')
		for i in range(3):
			self.print('\t\t'+str(opt_latt[i,:]))
		self.print('\t [optimize_lattice]: opt_vol='+str(opt_vol)+' (a0**3) **accepted** ( dV  ='+str(dV)+'(a0**3)	)')
    	#
		self.print('**** FINISHED STRAINED LATTIC CALCULATOR ********')
		return opt.x, opt_vol, opt_latt
	#~~~~~~~~~~~~~~~~~




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





def test():
	#	TEST THE CLASS
	dist_latt		=	FeMn_latt_OPT(1.1)
	tol_arr			=	np.array([1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-7,1e-15])
	
	print(r'  tol 	| x ($a_0$)	| vol ($a_0^3$) | $\deltaV (a_0^3)$	')
	print('----------------------------------')
	x_opt_lst	=	[]
	opt_vol_lst	=	[]
	opt_latt_lst=	[]
	dV_lst		=	[]


	for tol in tol_arr:
		x_opt, opt_vol, opt_latt	=	dist_latt.optimize_lattice(tol)
		dV 				=	np.abs(	dist_latt.vol0	- opt_vol)
		#
		x_opt_lst.append(		x_opt		)
		opt_vol_lst.append(		opt_vol		)
		opt_latt_lst.append(	opt_latt	)
		dV_lst.append(			dV 			)
		#
		print('{}{:e}{}{:f}{}{:f}{}{:f}'.format(' ',tol,' | ',x_opt_lst[-1]," | ",opt_vol_lst[-1]," | ",dV_lst[-1]))
	print("")
	print("opt latt:")
	print(opt_latt_lst[-1])


#		UNCOMMENT TO TEST CLASS
#test()





