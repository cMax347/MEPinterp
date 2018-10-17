import numpy as np
import re
import matplotlib.pyplot as plt

au_to_ev	= 27.21139


def get_high_symm_points(kpt_file):
	#high symmetry points have a descriptive character in their 5th column
	# these characters are read here
	#get the labels of each k-point
	kpt_label = []
	with open(kpt_file,'r',) as original_k:
		for idx, line in enumerate(original_k):
			if idx > 0:		#skip header
				line_label = ''
				if len(line[41:-1])>0:							#read everything after 4th column which ends at index 40 (hardcoded :ugly:)
					line_label = line[41:-1].replace(" ","")	#remove whitespace
				kpt_label.append(line_label)

	#find high symmertry points in labels 
	high_symm_points = []
	for idx,lbl in enumerate(kpt_label):
		if len(lbl) >0:
			my_lbl = lbl 
			if 'Gamma' in lbl:
				my_lbl = '$\Gamma$'
			high_symm_points.append([idx,my_lbl])
	return high_symm_points


def sort_energies(en_data):
	#transform the 1D en_data list into a 2D list
	#where en_plot[1:nkpts][1:nBands]
	en_plot = []
	k_old	= []
	for idx, en_val	in enumerate(en_data):
		kpt = en_val[1:4]
		en 	= en_val[4]

		#init k_old & append first value
		if idx == 0:
			k_old	= kpt
			en_plot.append( [])

		#in case of new kpt jump in first dimension
		if np.linalg.norm(kpt -k_old) > 1e-8:
			en_plot.append([])
			k_old = kpt
			#print('new k_old=',k_old,' idx=',idx)
		#append value to second dimension
		en_plot[-1].append(en)

	return en_plot




def plot_bandstruct(kpt_file, en_file, pdf_out_file, label_size=14, y_tick_size=12, plot_in_ev=False):

	#kpt_data = np.genfromtxt(kpt_file,skip_header=1,dtype=(float,float,float,float,str), missing_values='',filling_values='none')

	kpt_data 	= np.genfromtxt(kpt_file, skip_header=1, usecols= (0,1,2)	)

	en_data		= np.genfromtxt(en_file, skip_header=3, usecols=(0,1,2,3,4)	)


	print('found '+str(len(kpt_data))+' kpts')	
	high_symm_points	=	get_high_symm_points(kpt_file)	
	print('high symm pts:'+str(high_symm_points))

	#linspace for plotting
	k_plot		= np.linspace(	0.0,1.0,len(kpt_data)		)

	#map high symmetry points to the linspace
	k_labels	= []
	k_ticks		= []
	for symm_point in high_symm_points:
		k_ticks.append(		k_plot[		symm_point[0]	]		)
		k_labels.append(				symm_point[1]			)

	#sort energies by kpt and band index
	en_plot = sort_energies(en_data)

	print('raw en_plot:',en_plot)
	
	nBands =	len(en_plot[0])
	print('detected nBands='+str(nBands))



	#PLOTTING
	fig, ax  = plt.subplots(1,1) 

	#plot each band
	for band in range(nBands):
		en_band = []
		for idx,kpt in enumerate(k_plot):
			en_band.append(	en_plot[idx][band]	)
		if plot_in_ev:
			en_band	= np.array(en_band) * au_to_ev

		plt.plot(k_plot, en_band, '-',color='black')

	#x-axis
	ax.set_xlim([k_plot[0],k_plot[-1]])
	ax.set_xticks(k_ticks)
	ax.set_xticklabels(k_labels,fontsize=label_size)
	ax.grid(axis='x', alpha=.5, linewidth=.8,	color='black')


	#y-axis
	ax.set_yticks([-9,-6,-3,0,3,6])
	ax.set_ylim([-9.2,6.2])
	plt.tick_params(axis='y', which='major',left=True,right=True, direction='in',labelsize=y_tick_size)
	if plot_in_ev:
		plt.ylabel(r'$E \,(eV)$',fontsize=label_size)
	else:
		plt.ylabel(r'$E \,(E_h)$',fontsize=label_size)

	#save file
	plt.tight_layout()
	try:
		plt.savefig(pdf_out_file,bbox_inches='tight')
		print('saved band_structure: '+pdf_out_file)
	except:
		print('Error while saving the plot, try to show plot now in order to manually save it')
		plt.show()
	




def unit_test():
	plot_bandstruct('./kpts','./eBands.dat',
					'./bands.pdf',
					label_size=14,
					y_tick_size=12,
					plot_in_ev=False
				)





unit_test()








