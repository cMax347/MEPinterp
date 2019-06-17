import os
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
				stop
		else:
			print('[read_f90]: src folder"',self.src_dir,'" does not exist ')
			stop
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
		self.hw_lst				=	np.load(self.data_dir+'/hw_lst.npy')
		self.ef_lst				=	np.load(self.data_dir+'/ef_lst.npy')		
		self.occ_lst			=	np.load(self.data_dir+'/occ_lst.npy')
		self.smr_lst			=	np.load(self.data_dir+'/smr_lst.npy')
		#
		#	convert to eV
		self.hw_lst				=	(map_unit(	au_to_ev,	self.hw_lst		),	"eV")
		self.smr_lst			=	(map_unit(	au_to_ev,	self.smr_lst	),	"eV")
		self.ef_lst				=	(map_unit(	au_to_ev,	self.ef_lst		),	"eV")
		#


	def generic_read_response(self,fname):
		fpath	=	self.data_dir +'/'+fname
		#
		if os.path.isfile(fpath):
			try:
				tmp	=	np.load(fpath)
				return tmp
			except:
				print("[read_f90/generic_read_response]: could not read ",fpath)
		return None



	def read_2nd_photoC(self, SI=True):
		assumed_shape			=	(3,3,3,len(self.hw_lst[0]),len(self.smr_lst[0]),len(self.ef_lst[0]))
		#
		self.scndPhoto_data		=	np.load(self.data_dir+'/photoC_2nd.npy')	
		read_shape				=	self.scndPhoto_data.shape
		#
		#	TEST IF SHAPE IS CORRECT
		if (not same_shape(read_shape, assumed_shape)):
			print("[read_f90/read_2nd_photoC]: WARNING scndPhoto_data has wrong shape, ")
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if SI:
			#	^^^^^^^^^^^^^^^^^^^^^^
			#		CONVERT TO SI		[	e**2/hbar e/E_h	]_atomic	-> [A/V**2]
			#	
			# >>> [	e**2/hbar	]_atomic	=	2.434135×10^-4 	[	S	]_SI
			cond_quantum			=	(scpc_dic["atomic unit of charge"][0])**2
			cond_quantum			=	cond_quantum	/	(scpc_dic["atomic unit of action"][0])
			#
			# >>> [	e/E_h		]_atomic	=	0.03674932 		[	1/V	]_SI
			inv_field 				=	scpc_dic["atomic unit of charge"][0]/scpc_dic["Hartree energy"][0]
			#
			au_to_si			=	cond_quantum	*	inv_field		
			self.scndPhoto_data	=	map_unit(au_to_si, self.scndPhoto_data)
			#
			#
			self.scndPhoto_data	=	(self.scndPhoto_data,"A/V**2") 
			return self.scndPhoto_data
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#
		#
		#	ATOMIC UNITS
		self.scndPhoto_data	=	(self.scndPhoto_data, "e**2/hbar e/E_h")
		return self.scndPhoto_data
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		


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
#test_same_shape()



