import numpy as np
import datetime
import sys
import os
import matplotlib.pyplot as plt


from fortran_io			import 	read_real_tens_file
from fortran_io			import 	read_cmplx_tens_file

#class response_tensors:
#		def __init__(self):




class HW_probe:

	def __init__(self, root_dir):
		#derived attributes
		self.root_dir	= root_dir
		self.plot_dir	= self.root_dir+'/plots'
		self.n_bands	= 0

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
		
	

	def __del__(self):
		print("~")
		print('plotted hw probing, by\n\n')
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")


	

	def read_responses(self, out_dir, hw_idx):
		#
		hw_idx	=	hw_idx + 1	# Fortran index starts at 1
		id_str	=	'{:07d}'.format(hw_idx)
		print('read_responses by identyfing with id_str ="'+id_str+"'")
		#
		#	HALL LIKE
		ahc_tens		=	read_real_tens_file(	out_dir		+ 	'/ahc/ahc_tens.dat'			,			'ahc'		)		
		ahc_kubo_tens	=	read_cmplx_tens_file(	out_dir		+	'/ahc/ahc_velo.hw'	+ id_str,			'ahcVELO'	)
		ohc_kubo_tens	=	read_cmplx_tens_file(	out_dir		+	'/ahc/ohc_kubo.hw'	+ id_str,			'ohcVELO'	)
		#-------
		#
		#	OPTICAL
		SeCnd_opt_tens	=	read_real_tens_file(	 out_dir	+	'/opt/2nd_photo.hw'	+ id_str,			'2phC'		)
		optA_tens		=	read_cmplx_tens_file(	 out_dir	+	'/opt/opt_Asymm.hw'	+ id_str,			'optA'		)
		optS_tens		=	read_cmplx_tens_file(	 out_dir	+	'/opt/opt_Ssymm.hw'	+ id_str,			'optS'		)
		#-------
		#
		#	GYROTROPIC
		#-------
					#
		#
		#
		return 		ahc_tens, ahc_kubo_tens, ohc_kubo_tens, SeCnd_opt_tens, optA_tens, optS_tens


	def save_attach_subData(		self, work_dir, hw,
									ahc_tens, ahc_kubo_tens, ohc_kubo_tens,	
									SeCnd_opt_tens, optA_tens, optS_tens
							):
		#
		#only record if the container is not empty (i.e. the file was found and had good behaviour)
		#print("[save_attach_subData]: hw="+str(hw))
	
		#--------------------------
		#	HALL LIKE
		#
		if len(ahc_tens)	is 3:
			self.hw_ahc_data.append(				[hw, ahc_tens]			)
		else:
			print("[save_attach_subData]: 	WARNING  wrong ahc length in "+					str(work_dir)		)
			print("[save_attach_subData]: ahc_tens:"+str(ahc_tens))
		if len(ahc_kubo_tens) is 3:
			self.hw_ahc_kubo_data.append(			[hw, ahc_kubo_tens]		)
		else:
			print("[save_attach_subData]: 	WARNING  wrong ahc kubo length in "+			str(work_dir)		)
		if len(ohc_kubo_tens) is 3:
			self.hw_ohc_data.append(				[hw, ohc_kubo_tens]		)
		else:
			print("[save_attach_subData]: 	WARNING  wrong ohc kubo length in "+			str(work_dir)		)
		#--------------------------
		#	OPTICAL
		#
		if len(SeCnd_opt_tens)	is 3:
			self.hw_2ndPhoto_data.append(			[hw,	SeCnd_opt_tens]		)
		else:
			print("[save_attach_subData]: 	WARNING  wrong 2nd photo length in "+			str(work_dir)		)
		if len(optA_tens)	is 3:
			self.hw_optA_data.append(			[hw,	optA_tens]		)
		else:
			print("[save_attach_subData]: 	WARNING  wrong opt asymm tens length in "+		str(work_dir)		)
		if len(optS_tens)	is 3:
			self.hw_optS_data.append(			[hw,	optS_tens]		)
		else:
			print("[save_attach_subData]: 	WARNING  wrong opt ssymm tens length in "+		str(work_dir)		)





	def read_hw_lst(self,out_dir):
		self.hw_lst	=	np.genfromtxt(out_dir+'/hw_lst.txt', skip_header=1,usecols=(1))


	def read_ahc_lst(self,ahc_dir):
		self.hw_ahc_data.append(	)


	def read_hw(self):
		print("^")
		print("^")
		print("		TRAVERSE SUBDIRECOTRIES - COLLECT THE DATA	")
		print("-------------------------------------------------------------------------------")
		

		out_dir		=	self.root_dir+'/out'
		ahc_dir		=	out_dir+'/ahc'

		# READ HW LIST
		self.read_hw_lst(out_dir)


		#	loop hw_idx
		for hw_idx, hw_val in enumerate(self.hw_lst):
			print('[read_hw]: now try to read #hw ',hw_idx,"= ",hw_val," (eV)")
			ahc_tens, ahc_kubo_tens, ohc_kubo_tens, SeCnd_opt_tens, optA_tens, optS_tens	=	self.read_responses(out_dir,hw_idx)
			#
			self.save_attach_subData(	out_dir, hw_val, 	
										ahc_tens, ahc_kubo_tens, ohc_kubo_tens, 
										SeCnd_opt_tens, optA_tens, optS_tens
									)



		#	plot hall like neeeds
		# self.hw_lst:
		# self.hw_ahc_data:
		# self.hw_ahc_kubo_data:
		# self.hw_ohc_data:


		##	LOOP DIRECTORIES
		#for entry in os.scandir(self.root_dir):
		#	if entry.is_dir and 'plots' not in entry.path:
		#		#
		#		#
		#		work_dir	=	entry.path
		#		hw 		=	float(work_dir.split("hw")[1])	
		#		#
		#		#	READ DATA FROM NEW SUBDIR
		#		print("[read_hw]: found subdir ="+str(entry.path)+' assoc. hw='+str(hw))
		#		mep_tens, mep_cs, mep_lc, mep_ic,	ahc_tens, ahc_kubo_tens, ohc_kubo_tens,	SeCnd_opt_tens, optA_tens, optS_tens	=	self.read_subdirectory(work_dir)
		#		#
		#		self.save_attach_subData( work_dir,	hw,ahc_tens, ahc_kubo_tens, ohc_kubo_tens,	SeCnd_opt_tens, optA_tens, optS_tens	)
		#		#
		#		print("----------------------------------------------------------------------")
		#		print("")
		##
		##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		## ~					 ~
		##	SORT BY HW VALUE
		## ~					 ~
		##--
		#self.hw_lst				=	sorted( self.hw_lst						)
		
		##
		#self.hw_cs_data			=	sorted(	self.hw_cs_data					)
		#self.hw_lc_data			=	sorted(	self.hw_lc_data					)
		#self.hw_ic_data			=	sorted(	self.hw_ic_data					)
		##
		#self.hw_ahc_data		=	sorted(	self.hw_ahc_data				)
		#self.hw_ahc_kubo_data	=	sorted(	self.hw_ahc_kubo_data			)
		#self.hw_ohc_data		=	sorted(	self.hw_ohc_data				)
		##
		#self.hw_2ndPhoto_data	=	sorted(	self.hw_2ndPhoto_data			)
		#self.hw_optA_data		=	sorted(	self.hw_optA_data				)
		#self.hw_optS_data		=	sorted(	self.hw_optS_data				)
		##~~~~~~~~~~~~~~~~~~~~~~~~
		print("[read_hw]: hw list = ",self.hw_lst)








	def set_hall_units(self,units,scale):
		unit_dsc	=	"atomic units"
		unit_str	=	r'$e^2$/ ($\hbar a_0$)'
		#
		au_to_S 		=	3.874046 * 1e-5			#	(e**2/h)	-> (S) Siemens
		au_to_cm		=	5.2917721 * 1e-9


		au_to_S_cm		=	au_to_S	/ au_to_cm
		#
		if units == "scale":
			unit_dsc	=	"use the scale given as function argument (sort of a wildcard)"
			unit_str	=	"-"
		elif units == "SI":
			scale		=	au_to_S_cm
			unit_str		=	"[S/cm]"
			unit_dsc		=	"SI units"	
		elif units == "wx":		 	
			scale			=	au_to_S_cm	/ 100.
			unit_str		=	"[arb. u.]"
			unit_dsc		=	"Units used by wanxiang in his paper. this should be the SI value divided by 100"	
		#
		print('[set_hall_units]:  chooen units "'+units+'" with dim '+unit_str+'" and  descriptor: "'+unit_dsc+'" '	)
		#
		return scale, unit_str, unit_dsc 



			

	def plot_hall_like(		self, units='au', scale=1.0, 
							plot_ahc=True, plot_ahc_kubo= True, plot_ohc=True, 
							label_size=14, xtick_size=12, ytick_size=12,
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

		scale, unit_str, unit_dsc	=	self.set_hall_units(units,scale)

		for hw in self.hw_lst:
			hw_plot.append(hw)
		for ahc_tens in self.hw_ahc_data:
			ahc_data.append(						ahc_tens[1]				)
		for ahc_kubo_tens in self.hw_ahc_kubo_data:
			ahc_kubo_data.append(					ahc_kubo_tens[1]		)
		for ohc_tens in self.hw_ohc_data:
			ohc_kubo_data.append(					ohc_tens[1]				)
		#
		#
		
		#
		#
		dim_str	= []
		dim_str.append('x')
		dim_str.append('y')
		dim_str.append('z')

		#	uncomment for debugging purposes
		#print("[plot_hall_like]: hw list: 	", self.hw_tot_data	)
		#print("[plot_hall_like]: ahc list: 	", self.hw_ahc_data)
		#print("[plot_hall_like]: ahc kubo  list: 	", self.hw_ahc_kubo_data)
		#print("[plot_hall_like]: ohc kubo  list: 	", self.hw_ohc_data)
		hw_max	=	max(self.hw_lst)
		hw_min	=	min(self.hw_lst)


		#	color code for the AHC plot
		ahc_col		=	'black'
		ohc_wx_col	=	'blue'
		ohc_w90_col	=	'orange'	

		#
		#plot all 9 components of mep tensor
		for i in range(0,3):
			for j in range(0,3):
				#COLLECT a_ij(phi=0:n_phi)
				zero_plot			= []
				ahc_plot		 	= []
				RE_ahc_kubo_plot	= []
				IM_ahc_kubo_plot	= []
				RE_ohc_kubo_plot	= []
				IM_ohc_kubo_plot	= []
				#				
				for ahc_tens in ahc_data:
					ahc_plot.append(						scale * np.real(	ahc_tens[i][j]		)				)
					zero_plot.append(								0												)
				for ahc_K_tens in ahc_kubo_data:
					RE_ahc_kubo_plot.append(				scale * np.real(	ahc_K_tens[i][j]	)				)
					IM_ahc_kubo_plot.append(				scale * np.imag(	ahc_K_tens[i][j]	)				)
				for ohc_tens in ohc_kubo_data:
					RE_ohc_kubo_plot.append(				scale * np.real(	ohc_tens[i][j]		)				)
					IM_ohc_kubo_plot.append(				scale * np.imag(	ohc_tens[i][j]		)				)
				

				#
				#
				#do PLOT
				fig, ax  = plt.subplots(2,1, sharex=True)
				ax[0].set_title('optical Hall conductivity')
				#
				#		REAL PART
				#
				if plot_ahc:
					try:
						ax[0].plot(hw_plot, ahc_plot,'-', color=ahc_col,label=r'$ \: \: \:\: \;  \sigma ^{\mathrm{AHC}}$')
					except:
						print("[plot_hall_like]: 	WARNING could not plot AHC tensor")
				if plot_ahc_kubo:
					try:
						ax[0].plot(hw_plot, RE_ahc_kubo_plot,		'v-', 	color=ohc_wx_col, markersize=marker_size,		label=r'Wx'		)
					except:
						print("[plot_hall_like]: 	WARNING could not plot RE AHC_Kubo (wann guide) tensor")
				if plot_ohc:
					try:
						ax[0].plot(hw_plot, RE_ohc_kubo_plot,		'*-', 	color=ohc_w90_col, markersize=marker_size,		label=r'w90'		)
					except:
						print("[plot_hall_like]: 	WARNING could not plot RE AHC_Kubo (wanxiang) tensor")
				#
				#
				#		IMAGINARY PART
				#
				if plot_ahc_kubo:
					try:
						ax[1].plot(hw_plot, zero_plot,'-',color='black')
						ax[1].plot(hw_plot, IM_ahc_kubo_plot,		'v-', 	color=ohc_wx_col, markersize=marker_size,		label=r'Wx'		)
					except:
						print("[plot_hall_like]: 	WARNING could not plot IM AHC_Kubo (wann guide) tensor")
				if plot_ohc:
					try:
						ax[1].plot(hw_plot, IM_ohc_kubo_plot,		'*-', 	color=ohc_w90_col, markersize=marker_size,		label=r'w90'		)
					except:
						print("[plot_hall_like]: 	WARNING could not plot IM OHC_Kubo (wanxiang) tensor")
				#
				#
				#	LABELS & TITLES
				#
				ax[0].set_xticks(np.arange(hw_min, hw_max+1.0, 0.1), minor=True)
				try:	
					ax[0].set_xlim([hw_min, hw_max])
					#ax[0].set_ylim([- re_bound , re_bound  ])
					#ax[1].set_ylim([- im_bound , im_bound  ])
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
				ax[0].tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)
				ax[0].tick_params(axis='x',which='both', direction='in',labelsize=xtick_size)
				ax[1].tick_params(axis='x',which='both', direction='in',labelsize=xtick_size, top=True)
				ax[1].tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)
				#
				plt.legend()
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






def plot_hw(root_dir):
	#	use the above class in here to plot data in folder root_dir
	print('[plot_hw]:	hello there')
	if os.path.isdir(root_dir):
		myTest	= HW_probe(root_dir)
		print('[plot_hw]:	initalized plotting of '+str(root_dir))
		#
		myTest.read_hw()
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
									plot_ahc		=		True		, 
									plot_ahc_kubo	= 		True		, 
									plot_ohc		=		True		, 
									label_size=14, xtick_size=12, ytick_size=12, marker_size=2,
									re_bound		=	35,
									im_bound		=	10
							)
		print("...")
		print('[plot_hw]:	plotted Hall like tensors')
		#	~~~~~~~~~~~~~~~~~~~~~~~~
		#
		myTest.plot_opt()
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
	plot_hw(root_dir)



















###	OLD STUFF (DEPRECATED, MEP IS NOT FREQUENCY DEPENDENT, MAYBE THE PLOTTING ROUTINE COULD STILL BE REUSED FOR OTHER PLOTS SO KEEP IT FOR NOW)

#	def set_mep_units(self, units,scale):
#		unit_conv	=	1.0
#		unit_str	=	r'[$e^2$/$\hbar$]'
#		unit_dsc	=	'raw atomic units'
#		#
#		if units 	== 				'SI':
#			scale		=	2.434135e-4
#			unit_str	=	'[S]'
#			unit_dsc	=	'MEP in SI has units of conductivity [S]=[A/V]'
#		elif units	==			'cgs':
#			scale		=	0.09170124
#			unit_str	=	'[-]'
#			unit_dsc	=	'MEP is dimensionless quantity in cgs unit system'	
#		#
#		print("[set_mep_units]: MEP	unit id '"+units+"' was interpreted as "+unit_str+" the scale factor will be "+str(unit_conv))
#		#
#		return scale, unit_str, unit_dsc 	 
#
#	def plot_mep(self, units='au', scale=1.0, plot_contributions=True,plot_band_res=False, label_size=14, xtick_size=12, ytick_size=12):
#		print("^")
#		print("^")
#		print("-------------------------------------------------------------------------------")	
#		print("		PLOT MEP")
#		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
#		print("[plot_mep]:	try to create folder where plots should go")
#		#
#		if not os.path.isdir(self.plot_dir):
#			try:
#				os.mkdir(self.plot_dir)
#			except OSError:
#				print('[plot_mep]:	Could not make directory ',self.plot_dir, '	(proably exists already)')
#			finally:
#				print("~")
#		else:
#			print('[plot_mep]: ' + self.plot_dir + "	exists already! (WARNING older plots might be overwriten)")
#		#
#		hw_plot_bands	= []
#		hw_plot 		= []
#		#
#		mep_tot_data 	= []
#		mep_cs_data		= []
#		mep_lc_data		= []
#		mep_ic_data		= []
#		mep_band_data	= []
#		#
#		for hw in self.hw_lst:
#			hw_plot.append(hw)
#		print("[plot_mep]: hw_plot = ",hw_plot)
#
#		for data in self.hw_tot_data:
#			mep_tot_data.append(data[1])
#		#
#		if plot_band_res:
#			for n_hw,mep_band in enumerate(self.hw_bands_data):
#				#print('current PHI: #',n_phi)
#				phi_plot_bands.append(mep_band[0])
#				mep_band_data.append(mep_band[1])
#		#
#		scale, unit_str, unit_dsc 	=	self.set_mep_units(units,scale)
#		
#		#		
#		#
#		if plot_contributions:
#			print("[plot_mep]:	will plot Chern-Simons (CS) and  local & itinerant contributions ( lc & ic)")
#			for data in self.hw_cs_data:
#				mep_cs_data.append(data[1])
#			for data in self.hw_lc_data:
#				mep_lc_data.append(data[1])
#			for data in self.hw_ic_data:
#				mep_ic_data.append(data[1])
#		#
#		mep_max	=	.0
#		mep_min	=	.0
#		if len(mep_tot_data) > 0:
#			mep_max		= np.amax(mep_tot_data)
#			mep_min 	= np.amin(mep_tot_data)
#			mep_delta 	= (mep_max - mep_min) / 100.0
#			mep_max		= mep_max + mep_delta
#			mep_min		= mep_min - mep_delta
#		#
#		#
#		#
#		dim_str	= []
#		dim_str.append('x')
#		dim_str.append('y')
#		dim_str.append('z')
#		#
#		hw_max	=	max(self.hw_lst)
#		hw_min	=	min(self.hw_lst)
#		#
#		#plot all 9 components of mep tensor
#		for i in range(0,3):
#			for j in range(0,3):
#				#COLLECT a_ij(phi=0:n_phi)
#				mep_tot_plot 	= []
#				mep_cs_plot		= []
#				mep_lc_plot		= []
#				mep_ic_plot		= []
#				mep_band_plot	= []
#				for mep_tens in mep_tot_data:
#					mep_tot_plot.append(	scale * mep_tens[i][j]	)
#				if plot_contributions:
#					for mep_cs in mep_cs_data:
#						mep_cs_plot.append(		scale * mep_cs[i][j]	)
#					for mep_lc in mep_lc_data:
#						mep_lc_plot.append(		scale * mep_lc[i][j]	)
#					for mep_ic in mep_ic_data:
#						mep_ic_plot.append(		scale * mep_ic[i][j]	)
#
#				#collect band resolved data	
#				if plot_band_res:
#					#print('i=',i,' j=',j, '	band contributions:')
#					for phi, data in enumerate(mep_band_data):
#						#print('mep_band_data: #phi=',phi)
#						mep_band_plot.append([])
#						#print('---')
#						self.nn_bands	= 0
#						for idx,mep_band in enumerate(data):
#							self.nn_bands	=	self.nn_bands	+ 1
#							#print('band #',idx,'	data:',mep_band)
#							mep_band_plot[-1].append(	scale * mep_band[i][j])
#							#print('phi=',phi,' #band',idx, ' mep_band: ',mep_band[i][j])
#					#print('------')
#
#					#print("#bands	=	"+str(self.n_bands))
#					#print("mep_band_plot:",mep_band_plot)
#					#print('----------')
#					mep_band_final	= []
#
#					#n_bands	= 2
#					for band in range(self.n_bands):
#						mep_band_final.append([])
#						for hw in hw_plot_bands:
#							mep_band_final[-1].append(hw)
#
#					#print('mep_band_final container:',	mep_band_final)
#					
#					for hw,mep_bands in enumerate(	mep_band_plot	):
#						#print("phi=",phi,	'mep:')
#						#print(mep_bands)
#						for band, mep in enumerate(	mep_bands):
#							mep_band_final[band][hw]	=	mep
#
#					#print('mep_band_final values:',	mep_band_final)
#
#
#
#				#do PLOT
#				fig, ax  = plt.subplots(1,1) 
#				try:
#					plt.plot(hw_plot, mep_tot_plot,'-+', color='black',label="tot")
#				except:
#					print("[plot_mep]: ERROR cloud not plot mep_tot")
#
#				if plot_contributions:
#					try:
#						plt.plot(hw_plot, mep_cs_plot,		'--', 	color='red',		label="CS")
#					except:
#						print("[plot_mep]: ERROR cloud not plot mep_CS")
#					try:
#						plt.plot(hw_plot, mep_lc_plot,		'o-', 	color='darkblue',	label="LC")
#					except:
#						print("[plot_mep]: ERROR cloud not plot mep_LC")
#					try:
#						plt.plot(hw_plot, mep_ic_plot,		'^-', 	color='lightblue',	label="IC")
#					except:
#						print("[plot_mep]: ERROR cloud not plot mep_IC")
#
#				if plot_band_res:
#					try:
#						for band in range(self.n_bands):
#							plt.plot(hw_plot,	mep_band_final[band],	'x-', label=" #band "+str(band)		)
#					except:
#						print("[plot_mep]: WARNING some band contributions could not be plotted")
#				
#				#X-AXIS
#				plt.xlabel(r'$\hbar \omega$ (eV)',	fontsize=label_size)
#				ax.set_xlim([hw_min,hw_max])
#				#ax.set_xticks(	np.array([0,1,2])	,	minor=False	)
#				#ax.set_xticks(	np.array([0.5,1.5])	,	minor=True	)
#
#				#ax.set_xticklabels(np.array([r'$0$',r'$\pi$',r'$2\pi$']))
#				plt.tick_params(axis='x',which='major', direction='in',labelsize=xtick_size)
#				plt.tick_params(axis='x',which='minor', direction='in', labelsize=xtick_size)
#
#				#Y-AXIS
#				plt.ylabel(r'$\alpha_{'+dim_str[i]+dim_str[j]+'}$ '+unit_str,	fontsize=label_size)
#				#ax.set_ylim([mep_min,mep_max])
#				plt.tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)
#
#				if plot_contributions:
#					plt.legend()
#
#				plt.tight_layout()
#				plt.savefig(self.plot_dir+'/mep_'+dim_str[i]+dim_str[j]+'.pdf')
#				plt.close()
#				print('[plot_mep]:	finished processing '+dim_str[i]+dim_str[j]+' tensor, plot saved to: '+self.plot_dir+'/mep_'+dim_str[i]+dim_str[j]+'.pdf')
#		print("-------------------------------------------------------------------------------")
#		print("")
#		print("")	
#












