import os
import sys
import numpy as np
import scipy.constants as scpc 
#import scipy.constants.physical_constants as scpc_dic

au_to_ev	=	scpc.physical_constants["Hartree energy in eV"][0]
scpc_dic	=	scpc.physical_constants
dim_str		=	['x','y','z']


def map_unit(unit, arr):
	for idx,elem in enumerate(arr):
		arr[idx]	=	unit * elem
	#
	return arr


def same_shape(s1,s2):
	#
	d1	= len(s1)
	d2	= len(s2)
	#
	if (d1==d2):
		for dim in range(0,d1):
			if not(s1[dim]==s2[dim]):
				return False
		return True
	else:
		return False
	return False


def save_load_np(file):
	data = None
	try:
		data = np.load(file)
	except FileNotFoundError:
		print("[read_f90/load_np]: ERROR ",os.path.abspath(file)," does not exist")
		sys.exit(1)
	return data



class read_f90:
	def __init__(self,src):
		self.src_dir			=	os.path.abspath(	src		)
		#
		#	CHECK IF FOLDER EXISTS
		if os.path.isdir(self.src_dir):
			self.src_dir		=	os.path.abspath(		self.src_dir		)	
			if os.path.isdir(self.src_dir+'/out'):
				self.data_dir	=	os.path.abspath(	self.src_dir+'/out'		)
			else:
				print("[read_f90]: ERROR no 'out' folder provided in ",self.src_dir)
				sys.exit(1)
		else:
			print('[read_f90]: src folder"',self.src_dir,'" does not exist ')
			sys.exit(1)
		#
		#
		print("^^^^")
		print("[read_f90]: init")
		print("[read_f90]: src  dir: ",self.src_dir)
		print("[read_f90]: data dir: ",self.data_dir)
		#
		#	PARASPACE ARRAYS
		self.hw_lst				=	[]
		self.smr_lst			=	[]
		self.ef_lst				=	[]
		#
		self.hw_lst				=	save_load_np(self.data_dir+'/hw_lst.npy')
		self.ef_lst				=	save_load_np(self.data_dir+'/ef_lst.npy')		
		self.occ_lst			=	save_load_np(self.data_dir+'/occ_lst.npy')
		self.smr_lst			=	save_load_np(self.data_dir+'/smr_lst.npy')
		#
		#	convert to eV
		self.hw_lst				=	(map_unit(	au_to_ev,	self.hw_lst		),	"eV")
		self.smr_lst			=	(map_unit(	au_to_ev,	self.smr_lst	),	"eV")
		self.ef_lst				=	(map_unit(	au_to_ev,	self.ef_lst		),	"eV")
		#
		print("[read_f90]: PARASPACE INFO:")
		print("[read_f90]:  hw_lst==( ",self.hw_lst[0].shape,	', "',self.hw_lst[1],	'")')
		print("[read_f90]: smr_lst==( ",self.smr_lst[0].shape,	', "',self.smr_lst[1],	'")')
		print("[read_f90]:  ef_lst==( ",self.ef_lst[0].shape,	', "',self.ef_lst[1],	'")')

	
	def generic_read_response(self,fname):
		# this is a fallback function which tries to read a unsupported file 
		# allows for quick access of new tensors
		fpath	=	self.data_dir +'/'+fname
		#
		if os.path.isfile(fpath):
			try:
				tmp	=	save_load_np(fpath)
				return tmp
			except:
				print("[read_f90/generic_read_response]: could not read ",fpath)
		return None



	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#	ToDo:
	def read_ac_hall(self,dim=3,SI=True):
		assumed_shape			=	(3,3,len(self.hw_lst[0]),len(self.smr_lst[0]),len(self.ef_lst[0]))
		#
		self.ahc_dc_tens		=	save_load_np(self.data_dir+'/ahc_DC_tens.npy')
		self.ahc_ac_tens		=	save_load_np(self.data_dir+'/ahc_AC_tens.npy')
		#
		if dim==2:
			unit_str				=	r'$e^2$/ $\hbar$'
		elif dim==3:
			unit_str				=	r'$e^2$/ $\hbar a_0$'
		else:
			print("[read_ac_hall]:	WARNING unsupported dim=",dim," requested!")
		#
		
		#
		if SI: 
			unit_str			=	'S/m'
			# >>> [	e**2/hbar	]_atomic	=	2.434135×10^-4 	[	S	]_SI
			cond_quantum			=	(scpc_dic["atomic unit of charge"][0])**2
			cond_quantum			=	cond_quantum	/	(scpc_dic["atomic unit of action"][0])
			# >>> [ 1/ a_0 ]_atomic 		= 1.(5*10^-11)[ m ]
			length					=	1./(scpc_dic["Bohr radius"][0])

			au_to_si				=	cond_quantum*length
			print('[read_ac_hall]: au_to_si='+str(au_to_si))
			self.ahc_dc_tens		=	map_unit(au_to_si, self.ahc_dc_tens)
			self.ahc_ac_tens		=	map_unit(au_to_si, self.ahc_ac_tens)
		#
		self.ahc_dc_tens		=	(self.ahc_dc_tens, unit_str)
		self.ahc_ac_tens		=	(self.ahc_ac_tens, unit_str)
		#
		print("[read_f90]:  ahc_DC_tens==( ",self.ahc_dc_tens[0].shape,	', "',self.ahc_dc_tens[1],	'")')
		print("[read_f90]:  ahc_AC_tens==( ",self.ahc_ac_tens[0].shape,	', "',self.ahc_ac_tens[1],	'")')
		print("~~~~")
		#
		return self.ahc_dc_tens, self.ahc_ac_tens
	#
	def read_mep(self,SI=True):
		return None
	#

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	def read_2nd_photoC(self,dim=3, SI=True):
		assumed_shape			=	(3,3,3,len(self.hw_lst[0]),len(self.smr_lst[0]),len(self.ef_lst[0]))
		#
		self.scndPhoto_data		=	save_load_np(self.data_dir+'/photoC_2nd.npy')	
		read_shape				=	self.scndPhoto_data.shape
		#
		#	TEST IF SHAPE IS CORRECT
		if (not same_shape(read_shape, assumed_shape)):
			print("[read_f90/read_2nd_photoC]: WARNING scndPhoto_data has wrong shape, ")
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if SI:
			#	^^^^^^^^^^^^^^^^^^^^^^
			#		CONVERT TO SI		[ 	  e**2/hbar e/E_h	]_atomic,3D	-> [A/V**2]_3D
			#		CONVERT TO SI		[ a_0 e**2/hbar e/E_h	]_atomic,2D	-> [A m/V**2]_2D
			#	
			# >>> [	e**2/hbar	]_atomic	=	2.434135×10^-4 	[	S	]_SI
			cond_quantum			=	(scpc_dic["atomic unit of charge"][0])**2
			cond_quantum			=	cond_quantum	/	(scpc_dic["atomic unit of action"][0])
			#
			# >>> [	e/E_h		]_atomic	=	0.03674932 		[	1/V	]_SI
			inv_field 				=	scpc_dic["atomic unit of charge"][0]/scpc_dic["Hartree energy"][0]
			#
			au_to_si			=	cond_quantum	*	inv_field
			if dim==2:
				au_to_si		=	au_to_si * scpc_dic["atomic unit of length"][0]
				self.scndPhoto_data	=	map_unit(au_to_si, self.scndPhoto_data)
				self.scndPhoto_data	=	(self.scndPhoto_data,'A/V^ 2') 
			elif dim==3:
				self.scndPhoto_data	=	map_unit(au_to_si, self.scndPhoto_data)
				self.scndPhoto_data	=	(self.scndPhoto_data,'A/V^2') 
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		else:
			#	ATOMIC UNITS
			if dim==2:
				self.scndPhoto_data	=	(self.scndPhoto_data,'e^3 / (\hbar E_h) ')
			elif dim==3:
				self.scndPhoto_data	=	(self.scndPhoto_data, 'a_0 e^3 / (\hbar  E_h)')
		#
		print("[read_f90]:  scndPhoto_data==( ",self.scndPhoto_data[0].shape,	', "',self.scndPhoto_data[1],	'")')
		print("~~~~")
		return self.scndPhoto_data
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	
	def read_rashba_cfg(self,cfg_file='./rashba.cfg'):
		nMag 		= 	np.zeros(3)
		aR			=	0
		Vex			=	0
		#
		try:
			with open(cfg_file,'r') as cfg:
				for row,string in enumerate(cfg):
					if "=" in string:
						string	=	string.split("=")[1]
						string	=	string.split("#")[0]
						string	=	string.strip()
						print('striped string: "',string,'"')
					#
					if row==1:
						aR 		= 	float(	string	)
					elif row==2:
						Vex = 		float(	string	)#float(	string.split("=")[1]	)
					elif row==3:
						string_arr 	= 	string.split(" ")
						for idx, string in enumerate(string_arr):
							if idx<3:
								nMag[idx]	=	float(string)
		except FileNotFoundError:
			print("[read_rashba_cfg]: ERROR file ./rashba.cfg does not exist")
		#
		#	setup descriptive string
		descriptor	=	r' $\alpha_R=$'+'{:3.1f}'.format(aR)+r' $\mathrm{eV} \AA, \; V_{\mathrm{ex}}= $'+'{:3.1f}'.format(Vex)+r' $\mathrm{eV}$,'
		if(np.abs(Vex)>1e-3):
			if(	vec_is_parallel(nMag,ex)):
				descriptor	=	descriptor	+	r'$\;(\hat{\mathbf{n}}_{\mathrm{Mag}} \parallel 	\hat{\mathbf{e}}_x)$'
			if(	vec_is_parallel(nMag,ey)):
				descriptor	=	descriptor	+	r'$\;(\hat{\mathbf{n}}_{\mathrm{Mag}} \parallel 	\hat{\mathbf{e}}_y)$'
			if(	vec_is_parallel(nMag,ez)):
				descriptor	=	descriptor	+	r'$\;(\hat{\mathbf{n}}_{\mathrm{Mag}} \parallel 	\hat{\mathbf{e}}_z)$'
		#
		return aR, Vex, nMag, descriptor


def test():
	my_test	=	read_f90('./')
	#
	print("hw_lst: ",my_test.hw_lst)
	print("ef_lst: ",my_test.ef_lst)
	print("smr_lst: ",my_test.smr_lst)
	#
	photoC	=	my_test.read_2nd_photoC()
#test()



def test_same_shape():
	s0 = 	np.zeros((2,4)).shape
	s1 =	np.ones((2,4)).shape
	print('[test_same_shape]: test #0 passed? ',  same_shape(s0,s1)		)
	#
	s2 = 	np.ones((4,2)).shape
	print('[test_same_shape]: test #1 passed? ',not same_shape(s1,s2)	)
	#
	s3 =	np.ones((2,3,5)).shape
	print('[test_same_shape]: test #1 passed? ',not same_shape(s2,s3)	)
	print(scpc_dic["atomic unit of action"])
#test_same_shape()



