import numpy as np
import datetime
import sys
import os
import matplotlib.pyplot as plt
from random import randint
import matplotlib.cm as cm
import scipy.constants as scpc 
from read_f90 import read_f90

au_to_ev	=	scpc.physical_constants["Hartree energy in eV"][0]
print("scpc au_to_ev=",au_to_ev)

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    #cmap_name = base.name + str(N)
    #return base.from_list(cmap_name, color_list, N)
    return color_list





class plotter:

	def __init__(self, root_dir,dir_id,SI=True):
		#derived attributes
		self.root_dir	= root_dir
		self.data_dir	= root_dir+'/out'
		self.plot_dir	= self.root_dir+'/plots'
		self.dir_id		= dir_id
		self.n_bands	= 0
		#
		#
		#	containers
		#self.hw_lst				=	[]
		#self.ef_lst				=	[]
		#self.scndPhoto_data		=	[]
		#self.smr_lst			=	[]
#
#		##	read data
#		#self.hw_lst	=	np.load(self.data_dir+'/hw_lst.npy')		
#		#self.occ_lst=	np.load(self.data_dir+'/occ_lst.npy')
#		#self.smr_lst=	np.load(self.data_dir+'/smr_lst.npy')
#		#self.ef_lst	=	self.occ_lst[0][:]
#		#
#		##
#		#for idx, hw in enumerate(self.hw_lst):
#		#	self.hw_lst[idx]	=	hw	*	au_to_ev
#		#for idx, smr in enumerate(self.smr_lst):
#		#	self.smr_lst[idx]	=	smr	*	au_to_ev
#		##
		#self.scndPhoto_data	=	np.load(self.data_dir+'/photoC_2nd.npy')	
		
		self.data			=	read_f90(self.root_dir)
		raw					=	self.data.read_2nd_photoC(SI=SI)
		self.scndPhoto_data	=	raw[0]
		self.unit_str		=	'['+raw[1]+']'
		
		#
		np_arr				=	np.array(	self.scndPhoto_data)
		raw_shape			=	np_arr.shape
		#
		

	def __del__(self):
		print("~")
		print('plotted hw probing, by\n\n')
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")

			

	def plot_hall_like(		self, scale=1.0, phi_laser=1.0,
							plot_ahc=True, plot_ahc_kubo= True, plot_ohc=True, 
							line_width=1, label_size=14, xtick_size=12, ytick_size=12,
							marker_size=12,
							upper_bound=1, lower_bound=0,
							plot_legend=False
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
		dim_str	= []
		dim_str.append('x')
		dim_str.append('y')
		dim_str.append('z')
		#
		#
		smr_idx	=	1
		#	color code for the AHC plot
		colors 	= discrete_cmap(len(self.data.ef_lst[0]),	'cool')
		#
		#LOOP SPACIAL COMPONENTS OF TENSOR (make individual plot for each)
		for x in range(0,3):
			for i in range(0,3):
				for j in range(0,3):		
					#COLLECT a_ij(phi=0:n_phi)
					#do PLOT
					fig, ax  = plt.subplots(1,1, sharex=True)
					#				
					#
					title = ''
					plt.title(title)


					#	plot curv for each fermi level
					#for ef_idx, ef_val in enumerate(self.ef_lst):
					re_plot	=	[]
					im_plot	=	[]
					ef_idx =0
					ef_val	=	self.data.ef_lst[0][ef_idx]
					#
					# collect data for current plot
					for hw_idx, hw_val in enumerate(self.data.hw_lst[0]):
						re_plot.append(		np.real(	scale	*	self.scndPhoto_data[x][i][j][hw_idx][smr_idx][ef_idx]		))
						im_plot.append(		np.imag(	scale	*	self.scndPhoto_data[x][i][j][hw_idx][smr_idx][ef_idx]		))
					#
					#	plot
					ax.plot(self.data.hw_lst[0], re_plot,'-', color='deepskyblue'	,label='Re {:+4.2f}'.format(ef_val))
					ax.plot(self.data.hw_lst[0], im_plot,'-', color='darkorange'	,label='Im {:+4.2f}'.format(ef_val))

					#
					#
					#
					#	save plot when fermi levels are added
					#
					hw_max	=	max(self.data.hw_lst[0])
					hw_min	=	min(self.data.hw_lst[0])
					ax.set_xticks(np.arange(hw_min, hw_max+1.0, 0.1), minor=True)
					try:	
						ax.set_xlim([hw_min, hw_max])
						#ax.set_ylim([lower_bound, upper_bound  ])
						#
						ax.set(ylabel=r'$\sigma^{'+dim_str[x]+'}_{'+dim_str[i]+dim_str[j]+'}\;$' +	self.unit_str)
						ax.yaxis.label.set_size(label_size)
						
					except:
						print("[plot_hall_like]: labeling of plot failed")
					#
					plt.xlabel(r'$ \hbar \omega $ (eV)', 	fontsize=label_size)
					#
					#
					#ax.set_ylim([mep_min,mep_max])
					#ax[0].tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)
					#ax[0].tick_params(axis='x',which='both', direction='in',labelsize=xtick_size)
					#ax[1].tick_params(axis='x',which='both', direction='in',labelsize=xtick_size, top=True)
					#ax[1].tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)
					#
					if plot_legend:
						plt.legend(loc='lower right',title=r'$\mathrm{E}_f$ (eV); ($\Gamma$='+'{:6.3f}'.format(self.data.smr_lst[0][smr_idx])+self.data.smr_lst[1]+')')
					plt.tight_layout()
					#
					outFile_path	= self.plot_dir+'/Jphoto^'+dim_str[x]+'_'+dim_str[i]+dim_str[j]+'.pdf'
					plt.savefig(outFile_path)
					plt.close()
					#
					print('[plot_hall_like]:	finished processing '+dim_str[x]+'_'+dim_str[i]+dim_str[j]+' tensor, plot saved to: '+outFile_path	)
		print("-------------------------------------------------------------------------------")
		print("")
		print("")	


	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#	end of class 	HW_probe
#
#**************************************************************************************************************************************
#**************************************************************************************************************************************
#--------------------------------------------------------------------------------------------------------------------------------------
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv






def plot_scnd_photo():
	#	use the above class in here to plot data in folder root_dir
	print('[plot_scnd_photo]:	hello there')

	dir_id	=	"theta"

	root_dir	=	'.'

	if os.path.isdir(root_dir):
		#	read data
		myTest	= plotter(root_dir,dir_id,SI=False)
		#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		#
		#	PLOT RESPONSES 
		#
		#	~~~~~~~~~~~~~~~~~~~~~~~~
		#
		myTest.plot_hall_like(		scale			=		1.0			, 
									phi_laser		=		1.0			,	# = 1j
									plot_ahc		=		False		, 
									plot_ahc_kubo	= 		True		, 
									plot_ohc		=		False		, 
									line_width=1.5,label_size=14, xtick_size=12, ytick_size=12, marker_size=1.1,
									upper_bound		=	500,
									lower_bound		=	0,
									plot_legend=True
							)
		print("...")
		print('[plot_scnd_photo]:	plotted Hall like tensors')
		#	~~~~~~~~~~~~~~~~~~~~~~~~
		#
		#----------------------------------------------------------------------------------------------------------------------
	else:
		print('[plot_scnd_photo]:	ERROR '+str(root_dir)+'	seems to be non existing. please specify valid folder')








plot_scnd_photo()











