import numpy as np
import os


#def read_mep_file(file):
#	mep_tens = []
#	if os.path.exists(file):
#		mep_file_path	= file
#		with open(mep_file_path, 'r') as mep_file:
#			start = -10
#			for idx,line in enumerate(mep_file):
#				if 'begin mep' in line:
#					start = idx
#				if idx > start  and idx <= start + 3:
#					mep_tens.append(np.fromstring( line, dtype=np.float, sep=' ' ) )
#				if idx == start + 4 and 'end mep' not in line:
#					print('error at the end of reading mep tensor')
#
#		mep_tens = np.array(mep_tens)
#		print('read file '+mep_file_path)
#		if mep_tens.size is not 9:
#			print('WARNING issues with dimensionalty of MEP tensor')
#			print('mep_tens interpretation: '+str(mep_tens))
#	else:
#		print('could not find file '+file)
#	return mep_tens


def read_real_tens_file(file,id):
	tens = []
	if os.path.exists(file):
		tens_file_path	= file
		with open(tens_file_path, 'r') as tens_file:
			start = -10
			for idx,line in enumerate(tens_file):
				if 'begin '+id in line:
					start = idx
				if idx > start  and idx <= start + 3:
					row	= np.fromstring( line, dtype=np.float, sep=' ' )
					if row.size is 3:
						tens.append(row)
					else:
						print("[read_real_tens_file]	WARNING: FOUND ROW WHICH DID NOT HAVE 3 ENTRIES")
				if idx == start + 4 and 'end '+id not in line:
					print(		'error at the end of reading tensor in file '+str(file)		)
		#
		tens = np.array(tens)
		print(' > read real file '+tens_file_path)
		if tens.size is not 9:
			print('WARNING issues with dimensionalty of '+str(id)+' tensor')
			print(str(id)+'_tens interpretation: '+str(tens))
	else:
		print(' ! could not find file '+file)
	return tens


def read_cmplx_tens_file(file,id):
	tens	=	[]
	#
	if os.path.exists(file):
		#
		with open(file,'r') as tens_file:
			start = -10
			for idx, line in enumerate(tens_file):
				if 'begin '+id in line:
					start = idx
				if idx > start and idx <= start + 3:
					row = np.fromstring( line, dtype=np.float, sep=' ')
					if row.size is 6:
						tens.append(	np.array(		[	(row[0]+1j*row[1]),	(row[2]+1j*row[3]),		(row[4]+1j*row[5])	]		)		)
					else:
						print("[read_cmplx_tens_file]	WARNING: FOUND ROW WHICH DID NOT HAVE 6 ENTRIES ")
				if idx == start + 4 and 'end '+id not in line:
					print(		'error at the end of reading tensor in file '+str(file)		)
		#
		tens	=	np.array( tens ) 
		print(' > read cmplx file '+file)
	else:
		print(' ! could not find file '+file)
	return tens





