import numpy as np
import datetime
import sys
import os
import matplotlib.pyplot as plt


#from fortran_io			import 	read_real_tens_file
#from fortran_io			import 	read_cmplx_tens_file

#class response_tensors:
#		def __init__(self):




class HW_probe:

	def __init__(self, root_dir,dir_id):
		#derived attributes
		self.root_dir	= root_dir
		self.plot_dir	= self.root_dir+'/plots'
		self.dir_id		= dir_id
		self.n_bands	= 0
		#
		#
		dir_list		= os.listdir(self.root_dir)
		self.subdirs  	= []
		self.para_vals	= []
		for elem in dir_list:
			path = self.root_dir+'/'+elem
			if os.path.isdir(path):
				if (self.dir_id in elem):
					self.subdirs.append(path)
					self.para_vals.append(	elem.split(dir_id) [1]	)
		#
		#
		#	hw list & energies
		self.hw_lst				=	[]
		#
		#
		#	Hall like
		self.hw_ohc_data		=	[]
		self.hw_ahc_kubo_data	=	[]
		self.hw_ahc_data		=	[]
		#
		#	optical
		self.hw_2ndPhoto_data	=	[]
		self.hw_optA_data		=	[]
		self.hw_optS_data		=	[]
		#
		#	gyrotropic
		print("^")
		print("^")
		print("^")
		print("^")
		print("^^^^^^^^^^^^^^^	PLOTTING SCRIPT - SOUZA TB	-	RESPONSE OVER HOPPING PHASE	 ^^^^^^^^^^^^^^^")
		print("-------------------------------------------------------------------------------")
		print("~")
		print("[init]: will search for data in folder: "	+	self.root_dir	)
		print("[init]: will output to folder: "				+ 	self.plot_dir	)
		print("[init]: found data folders:")
		for idx, sub in enumerate(self.subdirs):
			print('[init]:		',sub, " intepreted as ",dir_id,'=',self.para_vals[idx])
	

	def __del__(self):
		print("~")
		print('plotted hw probing, by\n\n')
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")


	

#	def read_responses(self, out_dir, hw_idx):
#		#
#		hw_idx	=	hw_idx + 1	# Fortran index starts at 1
#		id_str	=	'{:07d}'.format(hw_idx)
#		print('read_responses by identyfing with id_str ="'+id_str+"'")
#		#
#		#	HALL LIKE
#		ahc_tens		=	read_real_tens_file(	out_dir		+ 	'/ahc/ahc_tens.dat'			,			'ahc'		)		
#		ahc_kubo_tens	=	read_cmplx_tens_file(	out_dir		+	'/ahc/ahc_velo.hw'	+ id_str,			'ahcVELO'	)
#		ohc_kubo_tens	=	read_cmplx_tens_file(	out_dir		+	'/ahc/ohc_kubo.hw'	+ id_str,			'ohcVELO'	)
#		#-------
#		#
#		#	OPTICAL
#		SeCnd_opt_tens	=	read_real_tens_file(	 out_dir	+	'/opt/2nd_photo.hw'	+ id_str,			'2phC'		)
#		optA_tens		=	read_cmplx_tens_file(	 out_dir	+	'/opt/opt_Asymm.hw'	+ id_str,			'optA'		)
#		optS_tens		=	read_cmplx_tens_file(	 out_dir	+	'/opt/opt_Ssymm.hw'	+ id_str,			'optS'		)
#		#-------
#		#
#		#	GYROTROPIC
#		#-------
#					#
#		#
#		#
#		return 		ahc_tens, ahc_kubo_tens, ohc_kubo_tens, SeCnd_opt_tens, optA_tens, optS_tens
#

#	def save_attach_subData(		self, work_dir,para_val, hw,
#									ahc_tens, ahc_kubo_tens, ohc_kubo_tens,	
#									SeCnd_opt_tens, optA_tens, optS_tens
#							):
#		#
#		#only record if the container is not empty (i.e. the file was found and had good behaviour)
#		#print("[save_attach_subData]: hw="+str(hw))
#	
#		#--------------------------
#		#	HALL LIKE
#		#
#		#
#		if len(ahc_tens)	is 3:
#			self.hw_ahc_data[-1].append(				[hw, ahc_tens]			)
#		else:
#			print("[save_attach_subData]: 	WARNING  wrong ahc length in "+					str(work_dir)		)
#			print("[save_attach_subData]: ahc_tens:"+str(ahc_tens))
#		if len(ahc_kubo_tens) is 3:
#			self.hw_ahc_kubo_data[-1].append(			[para_val, hw, ahc_kubo_tens]		)
#		else:
#			print("[save_attach_subData]: 	WARNING  wrong ahc kubo length in "+			str(work_dir)		)
#		if len(ohc_kubo_tens) is 3:
#			self.hw_ohc_data[-1].append(				[para_val, hw, ohc_kubo_tens]		)
#		else:
#			print("[save_attach_subData]: 	WARNING  wrong ohc kubo length in "+			str(work_dir)		)
#		#--------------------------
#		#	OPTICAL
#		#
#		if len(SeCnd_opt_tens)	is 3:
#			self.hw_2ndPhoto_data[-1].append(			[para_val, hw,	SeCnd_opt_tens]		)
#		else:
#			print("[save_attach_subData]: 	WARNING  wrong 2nd photo length in "+			str(work_dir)		)
#		if len(optA_tens)	is 3:
#			self.hw_optA_data[-1].append(			[para_val, hw,	optA_tens]		)
#		else:
#			print("[save_attach_subData]: 	WARNING  wrong opt asymm tens length in "+		str(work_dir)		)
#		if len(optS_tens)	is 3:
#			self.hw_optS_data[-1].append(			[para_val, hw,	optS_tens]		)
#		else:
#			print("[save_attach_subData]: 	WARNING  wrong opt ssymm tens length in "+		str(work_dir)		)
#




	def read_hw_lst(self,out_dir):
		#this_hw_lst				=	np.genfromtxt(out_dir+'/hw_lst.txt', skip_header=1,usecols=(1))
		this_hw_lst				=	np.load(out_dir+'/hw_lst.npy')
		self.hw_lst.append(	this_hw_lst	)


	def read_ahc_lst(self,ahc_dir):
		self.hw_ahc_data.append(	)


	def read_data(self):
		print("^")
		print("^")
		print("		TRAVERSE SUBDIRECOTRIES - COLLECT THE DATA	")
		print("-------------------------------------------------------------------------------")
		
		#if self.dir_id in curr_dir:
		#	delta  		=	curr_dir.split(self.dir_id) [1]
		#	print('found new subdir with delta = ', delta)
		#	#
		#	out_dir		=	curr_dir+'/out'
		#	ahc_dir		=	out_dir+'/ahc'
		#	#
		#	# READ HW LIST
		#	self.read_hw_lst(out_dir)
		


		#	get data from all subdirs
		for idx, path in enumerate(self.subdirs):
			out_dir	=	path + '/out'
			#
			# add container for this data set
			self.hw_ahc_data.append(		[])
			self.hw_ahc_kubo_data.append(	[])
			self.hw_ohc_data.append(		[])
			self.hw_2ndPhoto_data.append(	[])
			self.hw_optA_data.append(		[])
			self.hw_optS_data.append(		[])
			#
			#	get hw lst of out_dir
			self.read_hw_lst(out_dir)
			#	loop hw_idx
			#
			#	make a plot for this hw_lst 
			#for hw_idx, hw_val in enumerate(self.hw_lst[-1]):
			#	swag = 1
			#	#print('[read_hw]: now try to read #hw ',hw_idx,"= ",hw_val," (eV)")
			#	ahc_tens, ahc_kubo_tens, ohc_kubo_tens, SeCnd_opt_tens, optA_tens, optS_tens	=	self.read_responses(out_dir,hw_idx)
			#	#
			#	self.save_attach_subData(	out_dir, self.para_vals[idx], hw_val, 	
			#								ahc_tens, ahc_kubo_tens, ohc_kubo_tens, 
			#								SeCnd_opt_tens, optA_tens, optS_tens
			#
			self.hw_ahc_data[-1]		=	np.load(out_dir+'/ahc_tens.npy')
			self.hw_ahc_kubo_data[-1]	=	np.load(out_dir+'/ahcVELO.npy')							
			self.hw_ohc_data[-1]		=	np.load(out_dir+'/ohcVELO.npy')
			#
			print('[read_data]: read data from ',out_dir, ' (associated ',self.dir_id,'=',self.para_vals[idx],')')

	
		zipped	= zip(		self.para_vals, 
						self.hw_ahc_data, self.hw_ahc_kubo_data, self.hw_ohc_data, 
						self.hw_2ndPhoto_data, 
						self.hw_optA_data, self.hw_optS_data
					)
		sort	=	sorted(zipped)
		self.para_vals, self.hw_ahc_data, self.hw_ahc_kubo_data, self.hw_ohc_data, self.hw_2ndPhoto_data, self.hw_optA_data, self.hw_optS_data= map(list,zip(*sort))

		##~~~~~~~~~~~~~~~~~~~~~~~~
		




	def set_hall_units(self,units,scale):
		unit_dsc	=	"atomic units"
		unit_str	=	r'$e^2$/ ($\hbar a_0$)'
		#
		au_to_S 		=	2.4341348e-4	#3.874046 * 1e-5			#	(e**2/h)	-> (S) Siemens
		au_to_cm		=	5.2917721e-9	# * 1e-9


		au_to_S_cm		=	au_to_S	/ au_to_cm
		#
		if units == "scale":
			unit_dsc	=	"use the scale given as function argument (sort of a wildcard)"
			unit_str	=	"-"
		elif units == "SI":
			scale			=	scale * au_to_S_cm
			unit_str		=	"[S/cm]"
			unit_dsc		=	"SI units"	
		elif units == "wx":		 	
			scale			=	scale *au_to_S_cm	/ 100.
			unit_str		=	r'[$10^2$ S/cm]'
			unit_dsc		=	"Units used by wanxiang in his paper. this should be the SI value divided by 100"	
		#
		print('[set_hall_units]:  chooen units "'+units+'" with dim '+unit_str+'" and  descriptor: "'+unit_dsc+'" '	)
		#
		return scale, unit_str, unit_dsc 





			

	def plot_hall_like(		self, units='au', scale=1.0, 
							plot_ahc=True, plot_ahc_kubo= True, plot_ohc=True, 
							line_width=1, label_size=14, xtick_size=12, ytick_size=12,
							marker_size=12,
							re_bound=1, im_bound=1
					):
		print("^")
		print("^")
		print("-------------------------------------------------------------------------------")	
		print("		PLOT HALL LIKE")
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print("[plot_hall_like]:	try to create folder where plots should go")
		#
		#
		if not os.path.isdir(self.plot_dir):
			try:
				os.mkdir(self.plot_dir)
			except OSError:
				print('[plot_hall_like]:	Could not make directory ',self.plot_dir, '	(proably exists already)')
			finally:
				print("~")
		else:
			print('[plot_hall_like]: '+self.plot_dir+"	exists already! (WARNING older plots might be overwriten)")
		#		
		#
		hw_plot 		= []
		ahc_data	 	= []
		ahc_kubo_data	= []
		ohc_kubo_data	= []
		#
		#print('ac_kubo_data shape:',np.array(self.hw_ahc_kubo_data).shape)





		para_floats	=	[]
		for para_str in self.para_vals:
			para_floats.append(float(para_str))

		#print(self.hw_ahc_kubo_data)

		scale, unit_str, unit_dsc	=	self.set_hall_units(units,scale)


		print(self.hw_ahc_kubo_data[1])
		print('para_vals:',self.para_vals)
		


		#for idx,para_val in enumerate(self.para_vals):
		#	hw_plot.append([])
		#	ahc_data.append([])
		#	ahc_kubo_data.append([])
		#	ohc_kubo_data.append([])
		#	#
		#	for hw in self.hw_lst[idx]:
		#		hw_plot[-1].append(hw)
		#	for ahc_tens in self.hw_ahc_data[idx]:
		#		ahc_data[-1].append(						ahc_tens[1]				)
		#	for ahc_kubo_tens in self.hw_ahc_kubo_data[idx]:
		#		ahc_kubo_data[-1].append(					ahc_kubo_tens[2]		)
		#	for ohc_tens in self.hw_ohc_data[idx]:
		#		ohc_kubo_data[-1].append(					ohc_tens[2]				)

		#print(ahc_kubo_data)
		#
		#
	
		#
		#
		dim_str	= []
		dim_str.append('x')
		dim_str.append('y')
		dim_str.append('z')

		##	uncomment for debugging purposes
		#print("[plot_hall_like]: hw list: 	", self.hw_tot_data	)
		#print("[plot_hall_like]: ahc list: 	", self.hw_ahc_data)
		#print("[plot_hall_like]: ahc kubo  list: 	", self.hw_ahc_kubo_data)
		#print("[plot_hall_like]: ohc kubo  list: 	", self.hw_ohc_data)
		#hw_max	=	max(self.hw_lst)
		#hw_min	=	min(self.hw_lst)


		#	color code for the AHC plot
		zero_col=	'black'
		neg_col	=	'orange'
		pos_col	=	'blue'

		#
		#plot all 9 components of mep tensor
		for i in range(0,3):
			for j in range(0,3):
				#COLLECT a_ij(phi=0:n_phi)
				#do PLOT
				fig, ax  = plt.subplots(2,1, sharex=True)
				

				title = ''
				if plot_ahc:
					title = title + ' AHC '
				if plot_ahc_kubo:
					title = title + ' OHC-Wx'
				if plot_ohc:
					title = title + ' OHC-w90'
				ax[0].set_title(title)
				#
				for idx, para_val in enumerate(self.para_vals):	
					print('[plot_hall_like]: start with (',dim_str[i],',',dim_str[j],')',self.dir_id,'=',para_val)				
					ahc_plot		 	= []
					RE_ahc_kubo_plot	= []
					IM_ahc_kubo_plot	= []
					RE_ohc_kubo_plot	= []
					IM_ohc_kubo_plot	= []
					#
					#print('test shape:')
					#test = np.array(	self.hw_ahc_kubo_data[0]	)
					#print(test.shape)
					#for hw_idx, hw_val in enumerate(self.hw_lst):
					#	print('#hw_idx=',hw_idx)
					#	print(test[hw_idx]) 
					#print('3x3?!',self.hw_ahc_kubo_data[0])
					#ahc_plot.append(	self.ahc)	
					for hw_idx,hw_val in enumerate(self.hw_lst[idx]):
						RE_ahc_kubo_plot.append(	scale *	np.real(		np.array(self.hw_ahc_kubo_data[idx])[i][j][hw_idx]				))
						IM_ahc_kubo_plot.append(	scale *	np.imag(		np.array(self.hw_ahc_kubo_data[idx])[i][j][hw_idx]				))


					#print('idx=',idx,'RE_ahc_kubo_plot')
					#for elem in RE_ahc_kubo_plot:
					#	print(elem)
					#for ahc_tens in ahc_data[idx]:
					#	ahc_plot.append(						scale * np.real(	ahc_tens[i][j]		)				)
					#for ahc_K_tens in ahc_kubo_data[idx]:
					#	RE_ahc_kubo_plot.append(				scale * np.real(	ahc_K_tens[i][j]	)				)
					#	IM_ahc_kubo_plot.append(				scale * np.imag(	ahc_K_tens[i][j]	)				)
					#for ohc_tens in ohc_kubo_data[idx]:
					#	RE_ohc_kubo_plot.append(				scale * np.real(	ohc_tens[i][j]		)				)
					#	IM_ohc_kubo_plot.append(				scale * np.imag(	ohc_tens[i][j]		)				)
					#
					hw_max	=	max(self.hw_lst[idx])
					hw_min	=	min(self.hw_lst[idx])
					#
					#
					#	determine color based on para_val
					para_val = float(para_val)
					if abs( para_val-1.0 ) < 0.01:
						color	=	zero_col 
					elif para_val < 1.0:
						color	=	neg_col
					else:
						color	=	pos_col
					
					

					#	set the midrange curves to be dotted
					style='-'
					if (para_val > min(para_floats)) and (para_val < max(para_floats)):
						if( abs(para_val-1.0) > 0.01):
							style='-.' 


					label =r'$\delta$ '+'{:.2f}'.format(float(para_val))
					#		REAL PART
					#
					#if plot_ahc:
					#	#print('hw_plot=',hw_plot)
					#	#print('ahc_plot=',ahc_plot)
					#	try:
					#		ax[0].plot(hw_plot[0], ahc_plot,'-', color=color,label=label)
					#	except:
					#		print("[plot_hall_like]: 	WARNING could not plot AHC tensor")
					if plot_ahc_kubo:
						try:
							ax[0].plot(	self.hw_lst[idx], RE_ahc_kubo_plot,		
													marker='v',
													linestyle=style,
													linewidth=line_width, 	
													color=color ,	 
													markersize=marker_size,		
													label=label		
										)
						except:
							print("[plot_hall_like]: 	WARNING could not plot RE AHC_Kubo (wann guide) tensor")
					#$if plot_ohc:
					#$	try:
					#$		ax[0].plot(hw_plot[0], RE_ohc_kubo_plot,		
					#$									marker='*', 	
					#$									linestyle=style,
					#$									linewidth=line_width, 	
					#$									color=color, 
					#$									markersize=marker_size,		
					#$									label=label		)
					#$	except:
					#$		print("[plot_hall_like]: 	WARNING could not plot RE AHC_Kubo (wanxiang) tensor")
					#
					#
					#		IMAGINARY PART
					#
					if plot_ahc_kubo:
						try:
							ax[1].plot(self.hw_lst[idx], IM_ahc_kubo_plot,		
													marker='v',
													linestyle=style,
													linewidth=line_width, 	 	
													color=color, 
													markersize=marker_size,		
													label=label	
										)
						except:
							print("[plot_hall_like]: 	WARNING could not plot IM AHC_Kubo (wann guide) tensor")
					#if plot_ohc:
					#	try:
					#		ax[1].plot(hw_plot[0], IM_ohc_kubo_plot,		
					#									marker='*',
					#									linestyle=style,
					#									linewidth=line_width, 	 	
					#									color=color, 
					#									markersize=marker_size,		
					#									label=label	)
					#	except:
					#		print("[plot_hall_like]: 	WARNING could not plot IM OHC_Kubo (wanxiang) tensor")
				# ~~~~~~~~~~~~~~~~~
				#
				#
				#	LABELS & TITLES
				#
				ax[0].set_xticks(np.arange(hw_min, hw_max+1.0, 0.1), minor=True)
				try:	
					ax[0].set_xlim([hw_min, hw_max])
					ax[0].set_ylim([- re_bound , re_bound  ])
					ax[1].set_ylim([- im_bound , im_bound  ])
					#
					ax[0].set(ylabel=r'$\sigma^\mathrm{R}_{'+dim_str[i]+dim_str[j]+'}\; (\omega)\;$' +	unit_str)
					ax[0].yaxis.label.set_size(label_size)
					ax[1].set(ylabel=r'$\sigma^\mathrm{I}_{'+dim_str[i]+dim_str[j]+'}\; (\omega)\;$' +	unit_str)
					ax[1].yaxis.label.set_size(label_size)
				except:
					print("[plot_hall_like]: labeling of plot failed")
				#
				plt.xlabel(r'$ \hbar \omega $ (eV)',	fontsize=label_size)
				#
				#
				#ax.set_ylim([mep_min,mep_max])
				#ax[0].tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)
				#ax[0].tick_params(axis='x',which='both', direction='in',labelsize=xtick_size)
				#ax[1].tick_params(axis='x',which='both', direction='in',labelsize=xtick_size, top=True)
				#ax[1].tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)
				#
				plt.legend(loc='lower right')
				plt.tight_layout()
				fig.subplots_adjust(hspace=0.01)
				#
				outFile_path	= self.plot_dir+'/hall_'+dim_str[i]+dim_str[j]+'.pdf'
				plt.savefig(outFile_path)
				plt.close()
				#
				print('[plot_hall_like]:	finished processing '+dim_str[i]+dim_str[j]+' tensor, plot saved to: '+outFile_path	)
		print("-------------------------------------------------------------------------------")
		print("")
		print("")	






	def plot_opt(self):
		print("^")
		print("^")
		print("-------------------------------------------------------------------------------")	
		print("		PLOT HALL LIKE")
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print("[plot_opt]: 	WARNING this function is not implemented yet (ToDo!)	")
		#		MAYBE JUST ADD THIS TO HALL LIKE, THEN ITS EASIER TO COMPARE ALL DATA!!!!!!
		#

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#	end of class 	HW_probe
#
#**************************************************************************************************************************************
#**************************************************************************************************************************************
#--------------------------------------------------------------------------------------------------------------------------------------
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv






def plot_delta(root_dir):
	#	use the above class in here to plot data in folder root_dir
	print('[plot_hw]:	hello there')

	dir_id	=	"delta"



	if os.path.isdir(root_dir):
		# get data from sub dirs
		#
		# loop curr dir
		#	split dir by dir_id
		#d_val =	curr_dir.split(dir_id)[1]
		#
		#myTest	= HW_probe(curr_dir)

		myTest	= HW_probe(root_dir,dir_id)
		print('[plot_hw]:	initalized plotting of '+str(root_dir))
		#

		myTest.read_data()
		print('[plot_hw]:	read folder, start plotting....')
		print("..")
		print("..")
		#myTest.print_results_container()
		#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		#
		#	PLOT RESPONSES 
		#
		#	~~~~~~~~~~~~~~~~~~~~~~~~
		#
		myTest.plot_hall_like(		units			=		'wx'		, 
									scale			=		1.0			, 
									plot_ahc		=		False		, 
									plot_ahc_kubo	= 		True		, 
									plot_ohc		=		False		, 
									line_width=1.5,label_size=14, xtick_size=12, ytick_size=12, marker_size=1.1,
									re_bound		=	2.5,
									im_bound		=	5
							)
		print("...")
		print('[plot_hw]:	plotted Hall like tensors')
		#	~~~~~~~~~~~~~~~~~~~~~~~~
		#
		#myTest.plot_opt()
		print("...")
		print('[plot_hw]:	plotted optical cond. tensors')
		#----------------------------------------------------------------------------------------------------------------------
	else:
		print('[plot_hw]:	ERROR '+str(root_dir)+'	seems to be non existing. please specify valid folder')







#
#	GET FOLDER TO PLOT FROM STANDARD INPUT
#
if len(sys.argv) <	2:
	print("please pass a folder to plot")
else:
	root_dir	=	sys.argv[1]
	plot_delta(root_dir)











