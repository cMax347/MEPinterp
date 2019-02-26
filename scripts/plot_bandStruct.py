import numpy as np
import re
import os
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
		id	= int(en_val[0])
		kpt = en_val[1:4]
		en 	= en_val[4]
		#print('[sort_energies]:	new kpt:'+str(kpt)," new en:"+str(en))

		#init k_old & append first value
		if idx == 0:
			id_old 	= id
			k_old	= kpt
			en_plot.append( [])

		#in case of new kpt jump in first dimension
		if np.linalg.norm(kpt -k_old) > 1e-8:
			en_plot.append([])
			k_old = kpt
			if id == id_old:
				print("[plot_bandstruct/sort_energies]: WARNING unexpected new kpt found")
			#print('new k_old=',k_old,' idx=',idx)
		#append value to second dimension
		en_plot[-1].append(en)

	return en_plot



def get_en_plot(en_data):
	#sort energies by kpt and band index
	en_plot = sort_energies(en_data)
	nBands =	len(en_plot[0])
	#plot each band

	#
	return nBands, en_plot


def read_data(target_dir):
	k_plot		=	[]
	en_data		=	[]
	k_ticks		=	[]
	k_labels	=	[]
	#
	kpt_file	=	target_dir	+	'/kpts'
	en_file		=	target_dir 	+	'/out/eBands.dat'
	#
	if os.path.isfile(	kpt_file):
		kpt_data 	= np.genfromtxt(kpt_file, skip_header=1, usecols= (0,1,2)	)
		#linspace for plotting
		k_plot		= np.linspace(	0.0,1.0,len(kpt_data)		)
		print('[plot_bandstruct/read_data]: found '+str(len(kpt_data))+' kpts')
		high_symm_points	=	get_high_symm_points(kpt_file)
		print('[plot_bandstruct/read_data]: high symm pts:'+str(high_symm_points))
		#
		#
		for symm_point in high_symm_points:
			k_ticks.append(		k_plot[		symm_point[0]	]		)
			k_labels.append(				symm_point[1]			)
		#
	else:
		print("[plot_bandstruct/read_data]: ERROR did not find kpt_file "+kpt_file)
		stop
	if os.path.isfile(	en_file):
		en_data		= np.genfromtxt(en_file,  skip_header=3, usecols=(0,1,2,3,4)	)
	else:
		print("[plot_bandstruct/pread_data]: ERROR did not find en_file "+en_file)
		stop
	return 	k_plot, k_ticks, k_labels, en_data



def check_length(length, lst, default):
		orig_length	=	len(lst)
		if not (orig_length == length):
			print('[plot_bandstruct/check_length]: lst length was not of size '+str(length)+'	(lst has actual size '+str(orig_length)+')')
			print('[plot_bandstruct/check_length]: lst will be overwritten with default value = '+str(default))
			lst	=	[]
			for i in range(length):
				lst.append(default)
		if not (len(lst) ==  orig_length):
			print('[plot_bandstruct/check_length]: new lst size:'+str(len(lst)))
		return lst

def plot_bandstruct(target_dir_lst, id_str,id_formula,line_style, plot_color, pdf_out_file, label_size=14, y_tick_size=12, plot_in_ev=False):

	#kpt_data = np.genfromtxt(kpt_file,skip_header=1,dtype=(float,float,float,float,str), missing_values='',filling_values='none')


	print("[plot_bandstruct]: hello there! will search for data id: "+id_str)

	#this should be a unique identifier
	id_lst		=	[]
	line_style	=	check_length(len(target_dir_lst),	line_style,			"-"	)
	plot_color	=	check_length(len(target_dir_lst),	plot_color,		"black"	)


	#PLOTTING
	fig, ax  = plt.subplots(1,1)


	for dir_idx, next_dir in enumerate(target_dir_lst):
		if os.path.isdir(next_dir):
			print("[plot_bandstruct]:	NEW FOLDER FOUND	",next_dir," dir idx"+str(dir_idx))
			#
			#	print info on next_dir
			id_label	=	''
			try:
				id_lst.append(	float(next_dir.split(id_str)[1])	)
				id_label	=	id_formula+'='+'{:01.2f}'.format(id_lst[-1])
				print("[plot_bandstruct]: intepreted as "+str(id_str)+"="+id_label)
			except:
				print("[plot_bandstruct]: could not id the folder "+next_dir)
			#
			#
			k_plot, k_ticks, k_labels, en_data		=	read_data(next_dir)
			#
			#plot_color	=	'black'
			#line_style	=	'-'
			#
			#
			nBands,	en_plot	=	get_en_plot(en_data)
			print('[plot_bandstruct]: detected nBands='+str(nBands))
			for band in range(nBands):
				en_band = []
				for idx,kpt in enumerate(k_plot):
					en_band.append(	en_plot[idx][band]	)
				if plot_in_ev:
					en_band	= np.array(en_band) * au_to_ev
				if band == 0:
					plt.plot(k_plot, en_band, line_style[dir_idx],color=plot_color[dir_idx], label=id_label)
				else:
					plt.plot(k_plot, en_band, line_style[dir_idx], color=plot_color[dir_idx])
		else:
			print("[plot_bandstruct]: WARNING expected folder ",next_dir," was not found!")

	#x-axis
	ax.set_xlim([k_plot[0],k_plot[-1]])
	ax.set_xticks(k_ticks)
	ax.set_xticklabels(k_labels,fontsize=label_size)
	ax.grid(axis='x', alpha=.5, linewidth=.8,	color='black')




	#y-axis
	ax.set_ylim([-14,6.3])
	ax.set_yticks([0.],minor=True)
	ax.set_yticks([-12,-9,-6,-3,0,3,6],minor=False)
	plt.tick_params(axis='y', which='major',left=True,right=True, direction='in',labelsize=y_tick_size)
	plt.tick_params(axis='y', which='minor',left=True,right=True, direction='in',labelsize=y_tick_size-2)
	ax.grid(which='minor',axis='y')

	if plot_in_ev:
		plt.ylabel(r'$E \,(eV)$',fontsize=label_size)
	else:
		plt.ylabel(r'$E \,(E_h)$',fontsize=label_size)

	plt.legend(loc=(.26,.05),framealpha=1, shadow=False)

	#save file
	plt.tight_layout()
	#try:
	plt.savefig(pdf_out_file,bbox_inches='tight')
	print('saved band_structure: '+pdf_out_file)
	#except:
	#	print('Error while saving the plot, try to show plot now in order to manually save it')
        #		plt.show()





def unit_test():
	# 	out file
	pdf_target	=	'./bands.pdf'

	#	id data
	target_dir_lst	=[	'delta0.900',
						'delta0.950',
						'delta1.000',
						'delta1.050',
						'delta1.100'
					]
	id_str			=	'delta'
	id_label		=	r'$\delta$'
	#
	# 	colors
	min_col		=	'orangered'
	bas_col		=	'black'
	max_col		=	'deepskyblue'
	colors		=[min_col,min_col, bas_col, max_col, max_col]
	#
	#	line style
	prime		=	'-'
	opt			=	':'
	line_style	=	[prime, opt, prime, opt, prime]
	#
	#	
	plot_bandstruct(	target_dir_lst, 
						id_str, id_label,
						line_style, colors ,
						pdf_target, 
						label_size=14, y_tick_size=12, 
						plot_in_ev=True
					)
	
	#plot_bandstruct(['./'],'./kpts','./eBands.dat',
	#				'./bands.pdf',
	#				label_size=14,
	#				y_tick_size=12,
	#				plot_in_ev=False
	#			)
unit_test()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







