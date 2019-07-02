import sys
import numpy as np

import matplotlib.pyplot as plt






def check_if_fname_given():
	cnt 		= 	False
	data_file	=	"dummy"
	if len(sys.argv) >=	2:
		data_file	=	sys.argv[1]
		cnt 		=	True
		#
	return cnt, data_file


def read_file(fname):
	#
	#	try to read data with numpy 
	success			=	False
	raw_data		=	np.array([])
	try:
		raw_data	=	np.loadtxt(fname, usecols=(0,1,2) )
		success		=	True
		print("[read_file]: success! read data from file "+fname)
	finally:
		print("[read_file]: raw_data size	:	"+str(raw_data.size))
	#
	#	prepare data for plotting
	band_idx	=	[]
	k_val		=	[]
	en_val		=	[]
	k_val.append([])
	en_val.append([])
	print("[read_file]: start raw_data interpretation:")
	print("...")
	for idx,data_pt in enumerate(raw_data):
		#print('generic k val='+str(data_pt[0]),"	#band"+str(int(data_pt[1])),"	energy="+str(data_pt[2])			) 
		if len(band_idx) is 0:
			band_idx.append(int(data_pt[1]))
			#k_val.append([])
			#en_val.append([])
			print("		...initialized band_idx #"+str(band_idx[-1])+" array")
		elif int(data_pt[1]) != band_idx[-1]:
			band_idx.append(int(data_pt[1]))
			print("		...found new band #"+str(band_idx[-1]))
			k_val.append([])
			en_val.append([])
			#
		else:
			k_val[-1].append(			data_pt[0]			)
			en_val[-1].append(			data_pt[2]			)

		#
		#print(en_val)
	print("		...finished.")	
	print("[read_file]: Got a total of "+str(band_idx[-1])+" bands")
	return success, band_idx, k_val, en_val




def get_3qmodel_paras(q3_file):
	print("get_3qmodel_paras:	try to open "+str(q3_file) )
	with open(q3_file,'r') as inp_file:
		for idx,line in enumerate(inp_file):
			if idx==1:
				print(line)
				t1, t2, lmbda 	=	np.fromstring(line,dtype=float, count=3, sep=" ")
			if idx==2:
				phiA, thetaA	=	np.fromstring(line,dtype=float,count=2)
			if idx==3:
				phiB, thetaB	=	np.fromstring(line,dtype=float,count=2)	
			if idx==4:
				phiC, thetaC	=	np.fromstring(line,dtype=float,count=2)	
			if idx==5:
				phiD, thetaD	=	np.fromstring(line,dtype=float,count=2)		

	return t1, t2, lmbda, phiA, thetaA,  phiB, thetaB, phiC, thetaC, phiD, thetaD	


def plot_bandstruct(
				band_idx, k_val, en_val,
				pdf_out_file,
				q3_inp_file="./inp_params_3q",
				y_tick_size=12,
				label_size=12,
				title_size=14
	):
	#
	#INIT
	fig, ax  = plt.subplots(1,1) 

	#plot each band
	for band in band_idx:
		print("#"+str(band)+": got	"+str(		len(en_val[int(band)-1])		)+" kpts")
		#
		plt.plot(k_val[int(band)-1], en_val[int(band)-1], '-',color='black')
		#
		print("plotted band #"+str(band))

	# get 3q paras
	t1, t2, lmbda, phiA, thetaA,  phiB, thetaB, phiC, thetaC, phiD, thetaD	= get_3qmodel_paras(q3_inp_file)


	#
	#AESTHETICS
	#
	#x-axis
	#print("min k_val="+str(k_val[0][0]))
	#print("max k_val="+str(k_val[0][-1]))
	ax.set_xlim(			[ 	k_val[0][0], k_val[0][-1]	 ]					)		
	#ax.set_xticks(k_ticks)
	ax.set_xticklabels([],fontsize=label_size)
	#ax.grid(axis='x', alpha=.5, linewidth=.8,	color='black')
	#
	#y-axis
	#ax.set_yticks(	[-12, -9,-6,-3,0,3,6]	,	minor=False	)
	ax.set_yticks(	[0	]					,	minor=True		)
	ax.grid(axis='y',which='minor', alpha=.5, linewidth=.8,	color='black')
	##ax.set_ylim([-9.2,6.2])
	plt.tick_params(axis='y', which='major',left=True,right=True, direction='in',labelsize=y_tick_size)
	plt.ylabel(r'$E \,(eV)$',fontsize=label_size)
	#
	#title
	plt.title(r'3q FeMn:    $t_1=$'+str(t1)+r',$\: t_2=$'+str(t2)+r',$\: \lambda$='+str(lmbda), fontsize=title_size)

	#save file
	plt.tight_layout()
	try:
		plt.savefig(pdf_out_file,bbox_inches='tight')
		print('saved band_structure: '+pdf_out_file)
	except:
		print('Error while saving the plot, try to show plot now in order to manually save it')
		plt.show()
	








def plot_jan_bands(		pdf_out_file="./jans_bands.pdf",
						q3_inp_file="./inp_params_3q",
						y_tick_size=12,
						title_size=14
				):
	fname_provided, data_file	=	check_if_fname_given()
	#
	if fname_provided:
		print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
		print("|~~~~~~~~~~~~~~~~ READ DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~|")
		print("|-------------------------------------------------------|")
		f_read,	band_idx, k_val, en_val	=	read_file(data_file)
		if f_read:
			print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
			print("|~~~~~~~~~~~~~~~~ PLOT BANDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~|")
			print("|-------------------------------------------------------|")
			plot_bandstruct(
								band_idx, k_val, en_val,
								pdf_out_file,
								q3_inp_file,
								y_tick_size,
								title_size
							)
		else:
			print("could not read file"+str(data_file))


	else:
		print("please provide proper file name")



	
plot_jan_bands(
				pdf_out_file	=	"the_bands_of_jan.pdf"		,
				q3_inp_file		=	"./inp_params_3q"			,	
				y_tick_size		=	12							,
				title_size		=	13	
				)


