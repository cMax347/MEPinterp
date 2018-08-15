import numpy as np
import os


def read_mep_file(file):
	mep_tens = []
	if os.path.exists(file):
		mep_file_path	= file
		with open(mep_file_path, 'r') as mep_file:
			start = -10
			for idx,line in enumerate(mep_file):
				if 'begin mep' in line:
					start = idx
				if idx > start  and idx <= start + 3:
					mep_tens.append(np.fromstring( line, dtype=np.float, sep=' ' ) )
				if idx == start + 4 and 'end mep' not in line:
					print('error at the end of reading mep tensor')

		mep_tens = np.array(mep_tens)
		print('read file '+mep_file_path)
		if mep_tens.size is not 9:
			print('WARNING issues with dimensionalty of MEP tensor')
			print('mep_tens interpretation: '+str(mep_tens))
	else:
		print('could not find file '+file)
	return mep_tens
