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

def read_hw_lst(out_dir,fname='/hw_lst.txt'):
	hw_lst	=	np.genfromtxt(out_dir+fname, skip_header=1,usecols=(1))
	return hw_lst



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




def read_real_3rd_rank_tens_file(file,id, verbose=False):
	tens = []
	tens2= []
	#
	if os.path.exists(file):
		tens_file_path	= file
		with open(tens_file_path, 'r') as tens_file:
			start 	=	- 100
			end		=	- 100
			for idx,line in enumerate(tens_file):
				if 'begin '+id in line:
					start =	idx+1
					end	  = start + 9 + 3 
					if verbose:
						print('[read_real_3rd_rank_tens_file]: raw values:') 
				if idx == end:
					if not 'end '+id in line:
						print("[read_real_3rd_rank_tens_file]:unexpected end of file")
				#
				if idx > start:
					raw	=	np.fromstring(line, dtype=np.float, sep=' ')
					#	X-COMP
					if  idx <= start + 3:
						tens2.append(	raw	)
						if verbose:
							print('x _',idx-start,':',raw)
					# 	y init
					elif idx == start+4:
						tens.append(tens2)
						tens2	=	[]
					#	Y-COMP
					elif idx> start +4	and idx <=	start +7:
						tens2.append(	raw	)
						if verbose:
							print('y  _',idx-(start+4),':',raw)
					#	z init
					elif idx == start+8:
						tens.append(tens2)
						tens2	=	[]
					#	Z-COMP					
					elif idx> start+8 and idx <= start+11:
						tens2.append(	raw	)
						if verbose:
							print('z _',idx-(start+8),':',raw)
			# 
			tens.append(tens2)
		tens = np.array(tens)
	else:
		print('[read_real_3rd_rank_tens_file]: ERROR ',file,' does not exist')	
	#
	#
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





