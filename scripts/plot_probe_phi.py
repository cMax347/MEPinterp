import numpy as np
import datetime
import sys
import os
import matplotlib.pyplot as plt


from fortran_io			import 	read_real_tens_file





class Phi_probe:

	def __init__(self, root_dir):
		#derived attributes
		self.root_dir	= root_dir
		self.plot_dir	= self.root_dir+'/plots'
		self.n_bands	= 0

		print("~")
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")
		print("~")
		print("^^^^^^^^^^^^^^^	PLOTTING SCRIPT - SOUZA TB	-	RESPONSE OVER HOPPING PHASE	 ^^^^^^^^^^^^^^^")
		print("-------------------------------------------------------------------------------")
		print("~")
		print("will search for data in folder: "	+	self.root_dir	)
		print("will output to folder: "				+ 	self.plot_dir	)
		print("~")
		print("~")
		print("~")
		print("~")
		print("~")
		print("~")
		print("~")
		print("~")
		print("~")
		print("~")
		print("~")

	

	def __del__(self):
		print("~")
		print('plotted Phi probing, by\n\n')
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")


	def print_results_container(self):
		print("~")
		print("~")
		print("~")
		print("		COLLECTED DATA SUMMARY	")
		print("-------------------------------------------------------------------------------")
		for idx, phi_container in self.phi_tot_data:
			print("#phi=",idx)
			print(phi_container)
			print("----------------------------------")
		print("~")
		print("~")
		print("~")



	def read_phi(self):
		print("~")
		print("~")
		print("~")
		print("		TRAVERSE SUBDIRECOTRIES - COLLECT THE DATA	")
		print("-------------------------------------------------------------------------------")
		#
		#
		self.phi_tot_data 	= []
		self.phi_bands_data	=	[]
		self.phi_cs_data = []
		self.phi_lc_data = []
		self.phi_ic_data = []
		#
		for entry in os.scandir(self.root_dir):
			if entry.is_dir and 'plots' not in entry.path:
				#
				#
				work_dir	=	entry.path
				print("check subdir "+work_dir)
				phi 		=	work_dir.split("phi")[1]	
				print("found subdir ="+str(entry.path)+' assoc. phi='+phi)
				#
				#
				mep_file_path	= work_dir+'/out/mep/mep_tens.dat'
				mep_tens		= read_real_tens_file(mep_file_path,			'mep')
				#CHERN-SIMONS
				mep_file_path	= work_dir+'/out/mep/mep_cs.dat'
				mep_cs			= read_real_tens_file(mep_file_path,			'mep')
				#LOCAL
				mep_file_path	= work_dir+'/out/mep/mep_lc.dat'
				mep_lc		= read_real_tens_file(mep_file_path,				'mep')
				#ITINERANT
				mep_file_path	= work_dir+'/out/mep/mep_ic.dat'
				mep_ic		= read_real_tens_file(mep_file_path,				'mep')
				#
				search 		= True
				band 		= 1
				mep_bands	= []
				while search:
					mep_file_path 	=	work_dir+'/out/mep/mep_band.'+"{:07d}".format(band)
					tmp 			=	read_real_tens_file(mep_file_path,		'mep')
					if len(tmp)>0:
						#print('#band=',band, '	mep=',tmp)
						mep_bands.append(	tmp	)
						band = band + 1
					else:
						search	=	False
				self.n_bands = band -1 
				print("detected n_bands=",self.n_bands)

				print("phi=",phi," mep_bands=",mep_bands)
				#
				#
				#only record if the container is not empty (i.e. the file was found and had good behaviour)
				if len(mep_tens) is 3:
					self.phi_tot_data.append(		[phi, mep_tens	]			)
				else:
					print("skip mep_tens for phi="+str(phi))
					print("len(mep_tens)="+str(len(mep_tens)))
				if len(mep_cs) is 3:
					self.phi_cs_data.append(		[phi, mep_cs	]			)
				if len(mep_lc) is 3:
					self.phi_lc_data.append(		[phi, mep_lc	]			)
				if len(mep_ic) is 3:
					self.phi_ic_data.append(		[phi, mep_ic	]			)
				if len(mep_bands) >0:
					self.phi_bands_data.append(			[phi, mep_bands	]			)
				print("----------------------------------------------------------------------")
				print("")

		#			
		self.phi_tot_data 	= 	sorted(self.phi_tot_data)
		self.phi_cs_data	=	sorted(self.phi_cs_data)
		self.phi_lc_data	=	sorted(self.phi_lc_data)
		self.phi_ic_data	=	sorted(self.phi_ic_data)
		self.phi_bands_data	=	sorted(self.phi_bands_data)






	

	def plot_mep_over_phi(self, units='au', scale=1.0, plot_contributions=True,plot_band_res=False, label_size=14, xtick_size=12, ytick_size=12):
		print("		PLOTS")
		print("-------------------------------------------------------------------------------")	

		print("~")
		print("try to create folder where plots should go")


		if not os.path.isdir(self.plot_dir):
			try:
				os.mkdir(self.plot_dir)
			except OSError:
				print('Could not make directory ',self.plot_dir, '	(proably exists already)')
			finally:
				print("~")
		else:
			print(self.plot_dir+"	exists already! (WARNING older plots might be overwriten)")
		
		phi_plot_bands	= []
		phi_plot 		= []
		mep_tot_data 	= []
		mep_cs_data		= []
		mep_lc_data		= []
		mep_ic_data		= []
		mep_band_data	= []
		for data in self.phi_tot_data:
			phi_plot.append(float(data[0]))
			mep_tot_data.append(data[1])

		self.n_phi	=	0
		if plot_band_res:
			for n_phi,mep_band in enumerate(self.phi_bands_data):
				#print('current PHI: #',n_phi)
				phi_plot_bands.append(mep_band[0])
				mep_band_data.append(mep_band[1])
				self.n_phi	=	n_phi

		print("found ",self.n_phi,"	different phi values")

		
		print("~")
		print("fortran output is expected to be in atomic units	(dimensionless)")

		if units 	== 				'SI':
			unit_conv	=	2.434135e-4
			unit_str	=	'[S]'
		elif units	==			'cgs':
			unit_conv	=	7.297352726e-3 
			unit_str	=	'[cgs units - not implemented properly]'
		else							:
			unit_conv	=	1.0
			unit_str	=	'[-]'


		print("unit id '"+units+"' was interpreted as "+unit_str+" the scale factor will be "+str(unit_conv))
		print("")
		print("~")
		print("~")
		print("~")

	
		if plot_contributions:
			print("will plot Chern-Simons (CS) and  local & itinerant contributions ( lc & ic)")
			for data in self.phi_cs_data:
				mep_cs_data.append(data[1])
			for data in self.phi_lc_data:
				mep_lc_data.append(data[1])
			for data in self.phi_ic_data:
				mep_ic_data.append(data[1])

		print("~")
		mep_max		= np.amax(mep_tot_data)
		mep_min 	= np.amin(mep_tot_data)
		mep_delta 	= (mep_max - mep_min) / 100.0
		mep_max		= mep_max + mep_delta
		mep_min		= mep_min - mep_delta




		dim_str	= []
		dim_str.append('x')
		dim_str.append('y')
		dim_str.append('z')

		#plot all 9 components of mep tensor
		for i in range(0,3):
			for j in range(0,3):
				#COLLECT a_ij(phi=0:n_phi)
				mep_tot_plot 	= []
				mep_cs_plot		= []
				mep_lc_plot		= []
				mep_ic_plot		= []
				mep_band_plot	= []
				for mep_tens in mep_tot_data:
					mep_tot_plot.append(	scale * mep_tens[i][j]	)
				if plot_contributions:
					for mep_cs in mep_cs_data:
						mep_cs_plot.append(		scale * mep_cs[i][j]	)
					for mep_lc in mep_lc_data:
						mep_lc_plot.append(		scale * mep_lc[i][j]	)
					for mep_ic in mep_ic_data:
						mep_ic_plot.append(		scale * mep_ic[i][j]	)



				#collect band resolved data	
				if plot_band_res:
					#print('i=',i,' j=',j, '	band contributions:')
					for phi, data in enumerate(mep_band_data):
						#print('mep_band_data: #phi=',phi)
						mep_band_plot.append([])
						#print('---')
						self.nn_bands	= 0
						for idx,mep_band in enumerate(data):
							self.nn_bands	=	self.nn_bands	+ 1
							#print('band #',idx,'	data:',mep_band)
							mep_band_plot[-1].append(	scale * mep_band[i][j])
							#print('phi=',phi,' #band',idx, ' mep_band: ',mep_band[i][j])
					#print('------')

					#print("#bands	=	"+str(self.n_bands))
					#print("mep_band_plot:",mep_band_plot)
					#print('----------')
					mep_band_final	= []

					#n_bands	= 2
					for band in range(self.n_bands):
						mep_band_final.append([])
						for phi in phi_plot_bands:
							mep_band_final[-1].append(phi)

					#print('mep_band_final container:',	mep_band_final)
					
					for phi,mep_bands in enumerate(	mep_band_plot	):
						#print("phi=",phi,	'mep:')
						#print(mep_bands)
						for band, mep in enumerate(	mep_bands):
							mep_band_final[band][phi]	=	mep

					#print('mep_band_final values:',	mep_band_final)



				#do PLOT
				fig, ax  = plt.subplots(1,1) 
				plt.plot(phi_plot, mep_tot_plot,'-+', color='black',label="tot")

				if plot_contributions:
					plt.plot(phi_plot, mep_cs_plot,		'--', 	color='red',		label="CS")
					plt.plot(phi_plot, mep_lc_plot,		'o-', 	color='darkblue',	label="LC")
					plt.plot(phi_plot, mep_ic_plot,		'^-', 	color='lightblue',	label="IC")

				if plot_band_res:

					for band in range(self.n_bands):
						plt.plot(phi_plot,	mep_band_final[band],	'x-', label=" #band "+str(band)		)
				
				#X-AXIS
				plt.xlabel(r'$\varphi_{\mathcal{Hopp}}$',	fontsize=label_size)
				ax.set_xlim([0,2])
				ax.set_xticks(	np.array([0,1,2])	,	minor=False	)
				ax.set_xticks(	np.array([0.5,1.5])	,	minor=True	)

				ax.set_xticklabels(np.array([r'$0$',r'$\pi$',r'$2\pi$']))
				plt.tick_params(axis='x',which='major', direction='in',labelsize=xtick_size)
				plt.tick_params(axis='x',which='minor', direction='in', labelsize=xtick_size)

				#Y-AXIS
				plt.ylabel(r'$\alpha_{'+dim_str[i]+dim_str[j]+'}$ '+unit_str,	fontsize=label_size)
				#ax.set_ylim([mep_min,mep_max])
				plt.tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)

				if plot_contributions:
					plt.legend()

				plt.tight_layout()
				plt.savefig(self.plot_dir+'/mep_'+dim_str[i]+dim_str[j]+'.pdf')
				plt.close()
				print('finished processing '+dim_str[i]+dim_str[j]+' tensor, plot saved to: '+self.plot_dir+'/mep_'+dim_str[i]+dim_str[j]+'.pdf')


				

	
#**************************************************************************************************************************************
#**************************************************************************************************************************************
#**************************************************************************************************************************************
#**************************************************************************************************************************************
#**************************************************************************************************************************************







def plot_data(root_dir):
	myTest	= Phi_probe(root_dir)
	#
	myTest.read_phi()
	myTest.print_results_container()
	
	myTest.plot_mep_over_phi(		scale=1.0,
									units="au",
									plot_contributions=False,
									plot_band_res=True,
									label_size=14, xtick_size=12, ytick_size=12)






if len(sys.argv) <	2:
	print("please pass a folder to plot")
else:
	root_dir	=	sys.argv[1]
	plot_data(root_dir)







