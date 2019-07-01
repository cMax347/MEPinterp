import numpy as np
import datetime
import sys
import os
import matplotlib.pyplot as plt
from read_f90 import read_f90
#~
#~
#~
dim_str	= []
dim_str.append('x')
dim_str.append('y')
dim_str.append('z')
#~
class harry_plotter:
	def __init__(self,SI=True):
		#derived attributes
		self.root_dir	= os.getcwd()
		self.plot_dir	= os.path.abspath(	self.root_dir+'/plots'		)
		#
		#
		print("^")
		print("^")
		print("^")
		print("^")
		print("^^^^^^^^^^^^^^^	AC HALL PLOTTING SCRIPT ^^^^^^^^^^^^^^^")
		print("-------------------------------------------------------------------------------")
		print("try to read data")
		self.data								=	read_f90(self.root_dir)
		self.ahc_dc_tens, self.ahc_ac_tens		=	self.data.read_ac_hall(SI=SI)
		
		print("~")
		print("[init]: will output to folder: "				+ 	os.path.dirname(os.path.abspath(self.plot_dir)))
		#
		#
		#
		#
	#~~~~~~~~
	#
	#
	def __del__(self):
		print("~")
		print('plotted hw probing, by\n\n')
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")
	#~~~~~~~~
	#
	#
	def plot_ac_hall(		self, scale=1.0, scale_str='',
							plot_ahc=True, plot_ahc_kubo= True, plot_ohc=True, 
							line_width=1, label_size=14, xtick_size=12, ytick_size=12,
							marker_size=12,
							re_bound=1, im_bound=1
					):
		print("^")
		print("^")
		print("-------------------------------------------------------------------------------")	
		print("		PLOT AC HALL VS FREQUENCY")
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		print("[plot_ac_hall]:	try to create folder where plots should go")
		#
		#
		if not os.path.isdir(self.plot_dir):
			try:
				os.mkdir(self.plot_dir)
			except OSError:
				print('[plot_ac_hall]:	Could not make directory ',self.plot_dir, '	(proably exists already)')
			finally:
				print("~")
		else:
			print('[plot_ac_hall]: '+self.plot_dir+"	exists already! (WARNING older plots might be overwriten)")
		#		
		#
		hw_plot 		= []
		ahc_data	 	= []
		ahc_kubo_data	= []
		ohc_kubo_data	= []
		#
		unit_str		=	'('+scale_str + self.ahc_ac_tens[1]+')'
		print("unit str:",unit_str)
		#
		#
		ef_idx	=	2
		smr_idx	=	0
		#
		re_color = 'deepskyblue'
		im_color = 'gold'
		#
		re_label = r''
		im_label = r''
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
				ahc_plot		 	= []
				RE_ahc_kubo_plot	= []
				IM_ahc_kubo_plot	= []
				#
				for hw_idx,hw_val in enumerate(self.data.hw_lst[0]):
					ahc_plot.append(	scale * self.ahc_dc_tens[0][i][j][ef_idx])
					RE_ahc_kubo_plot.append(	scale *	np.real(		self.ahc_ac_tens[0][i][j][hw_idx][smr_idx][ef_idx]		))
					IM_ahc_kubo_plot.append(	scale *	np.imag(		self.ahc_ac_tens[0][i][j][hw_idx][smr_idx][ef_idx]		))
				#
				hw_max	=	max(self.data.hw_lst[0])
				hw_min	=	min(self.data.hw_lst[0])
				#
				#	set the midrange curves to be dotted
				style='-'
				#	DC REAL PART
				ax[0].plot(	self.data.hw_lst[0], ahc_plot,marker='v',linestyle=style,linewidth=line_width, color='black' ,	 markersize=marker_size)
				#	AC REAL PART
				ax[0].plot(	self.data.hw_lst[0], RE_ahc_kubo_plot,		
										marker='v',
										linestyle=style,
										linewidth=line_width, 	
										color=re_color ,	 
										markersize=marker_size,		
										label=re_label	+r'@$E_f=$'+str(self.data.ef_lst[0][ef_idx])+''	
							)
				#
				#	AC IMAG PART
				ax[1].plot(	self.data.hw_lst[0], IM_ahc_kubo_plot,		
										marker='v',
										linestyle=style,
										linewidth=line_width, 	
										color=im_color ,	 
										markersize=marker_size,		
										label=im_label+r'@$E_f=$'+'{:4.2f}'.format(self.data.ef_lst[0][ef_idx])+''			
							)
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
					print("[plot_ac_hall]: labeling of plot failed")
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
				print('[plot_ac_hall]:	finished processing '+dim_str[i]+dim_str[j]+' tensor, plot saved to: '+outFile_path	)
		print("-------------------------------------------------------------------------------")
		print("")
		print("[plot_ac_hall]: selected smearing: ",self.data.smr_lst[0][smr_idx]," eV; fermi_level:",self.data.ef_lst[0][ef_idx]," eV")
		print("")		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#	end of class 	HW_probe
#
#**************************************************************************************************************************************
#**************************************************************************************************************************************
#--------------------------------------------------------------------------------------------------------------------------------------
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv






def plot_ac_hall_tens(SI=True):
	#	use the above class in here to plot data in folder root_dir
	print('[plot_ac_hall_tens]:	hello there')

	dir_id	=	"delta"
	#
	myTest	= harry_plotter(SI=SI)
	#
	print('[plot_ac_hall_tens]:	read folder, start plotting....')
	print("..")
	print("..")
	#myTest.print_results_container()
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	#
	#	PLOT RESPONSES 
	#
	#	~~~~~~~~~~~~~~~~~~~~~~~~
	#
	myTest.plot_ac_hall(		scale			=		1		, 
								scale_str		=		r'',
								plot_ahc		=		False		, 
								plot_ahc_kubo	= 		True		, 
								plot_ohc		=		False		, 
								line_width=1.5,label_size=12, xtick_size=12, ytick_size=12, marker_size=1.1,
								re_bound		=	10,
								im_bound		=	10
						)
	print("...")
	print('[plot_ac_hall_tens]:	plotted Hall like tensors')
	#	~~~~~~~~~~~~~~~~~~~~~~~~
	#
	#myTest.plot_opt()
	print("...")
	print('[plot_ac_hall_tens]:	plotted AC hall tensor')
	#----------------------------------------------------------------------------------------------------------------------
#
plot_ac_hall_tens(SI=False)
#~~~~~~~~~~~~~











