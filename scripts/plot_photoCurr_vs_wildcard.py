import numpy as np
import datetime
import sys
import os
import matplotlib.pyplot as plt
from random import randint
import matplotlib.cm as cm
from read_f90 import read_f90
from laser import laser
#
dim_str	= []
dim_str.append('x')
dim_str.append('y')
dim_str.append('z')
#
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
#~~~~~~~~~~~~~~~~~~~




class merry_plotter:
	def __init__(self, root_dir,SI=True):
		#derived attributes
		self.root_dir	= root_dir
		self.data_dir	= root_dir+'/out'
		self.plot_dir	= self.root_dir+'/plots'
		self.n_bands	= 0
		self.SI_switch	= SI
		#
		#
		print("[merry_plotter]:	try to create folder where plots should go")
		if not os.path.isdir(self.plot_dir):
			try:
				os.mkdir(self.plot_dir)
			except OSError:
				print('[merry_plotter]:	Could not make directory ',self.plot_dir, '	(proably exists already)')
			finally:
				print("~")
		else:
			print('[merry_plotter]: '+self.plot_dir+"	exists already! (WARNING older plots might be overwriten)")
	#~
	#~
	def __del__(self):
		print("~")
		print('plotted hw probing, by\n\n')
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")
		print("-------------------------------------------------------------------------------")
	#~
	#~		
	def plot_photoC(		self, scale=1.0,scale_str='', 
							dim=3,
							x_axis_para_id=0,
							smr_idx=0, ef_idx=0,hw_idx=0,
							laser=None, laser_dir=2	,laser_lmbda=0,		 
							line_width=1, label_size=14, xtick_size=12, ytick_size=12,
							marker_size=12,
							upper_bound=1, lower_bound=0,
							plot_legend=False,
							title=''
					):
		#
		#	read photo c tens
		self.data			=	read_f90(self.root_dir)
		raw					=	self.data.read_2nd_photoC(dim=dim, SI=self.SI_switch)
		self.scndPhoto_data	=	raw[0]
		self.raw_unit		=	raw[1]
		#
		print("^")
		print("^")
		print("-------------------------------------------------------------------------------")	
		print("		PLOT 2nd order PHOTOCURRENT")
		print("              ~~~ VS ~~~  ")
		#
		if x_axis_para_id==0:
			para_str = 'smr'				
		elif x_axis_para_id==1:
			para_str = 'ef'	
		else:
			para_str = 'hw'	
		#
		print("                "+para_str.upper())
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		#
		self.unit_str	=		r'($'+scale_str+' '+self.raw_unit+laser.E0_units+laser.E0_units+r'$)'
		#
		if dim==3:
			new_unit_str	=		r'$(A / m^2)_{\mathrm{3D}}$'
		elif dim==2:
			new_unit_str	=		r'$(A / m)_{\mathrm{2D}}$'

		print("[plot_photoC]: unit info")
		print("")
		print("[plot_photoC]: raw unit str: ",self.unit_str)
		print("[plot_photoC]: WARNING interpreted as ",new_unit_str)
		print("")
		self.unit_str	=	new_unit_str
		#
		#	make sure smr_ and ef_idx are valid
		if smr_idx > len(self.data.smr_lst[0])-1:
			smr_idx = len(self.data.smr_lst[0])-1
			print("[plot_photoC]: WARNING smr_idx out of bounds was set, smr_idx set to ",smr_idx)
		if ef_idx > len(self.data.ef_lst[0])-1:
			ef_idx = len(self.data.ef_lst[0])-1
			print("[plot_photoC]: WARNING ef_idx out of bounds was set, ef_idx set to ",ef_idx)
		
		if hw_idx > len(self.data.hw_lst[0])-1:
			hw_idx = len(self.data.hw_lst[0])-1
			print("[plot_photoC]: WARNING hw_idx out of bounds was set, hw_idx set to ",hw_idx)
		#
		smr_val	=	self.data.smr_lst[0][smr_idx]
		ef_val	=	self.data.ef_lst[0][ef_idx]
		hw_val  = 	self.data.hw_lst[0][hw_idx]
		#
		#	color code for the AHC plot
		colors 	= discrete_cmap(len(self.data.ef_lst[0]),	'cool')
		#
		print("^")
		print("^")
		print("--------------------")	
		print("		LASER WARMUP")
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		laser.set_pol(laser_dir, laser_lmbda)
		print("\tlaser propagation dirct: ",dim_str[laser_dir])
		print("\tlaser pol: ",laser.pol)
		laser_E0	=	laser.get_field_strength()
		print("laser field strength E_0**2=",laser_E0**2,'	(V/m)**2')
		print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
		#
		dir_color	=	["deepskyblue","darkorange","purple"]
			# todo colorcode the current direcrtion

		#
		#LOOP SPACIAL COMPONENTS OF TENSOR (make individual plot for each)
		fig, ax  = plt.subplots(1,1, sharex=True)
		for x in range(0,dim):
			re_plot	=	[]
			#
			# collect data for current plot
			#^^^^^^^^^^^^^^^^^^^^^^
			#^^^^^^^^^^^^^^^^^^^^^^
			#^^^^^^^^^^^^^^^^^^^^^^
			#^^^^^^^^^^^^^^^^^^^^^^
			if x_axis_para_id==0:		#	VS SMEARING
				#collect
				for smr_idx, smr_val in enumerate(self.data.smr_lst[0]):
					re_plot.append(0)
					for i in range(0,dim):
						for j in range(0,dim):
							phi_laser	=	laser.get_phase(i,j)
							re_plot[-1]		= 		re_plot[-1] + np.real(	scale *	phi_laser *	self.scndPhoto_data[x][i][j][hw_idx][smr_idx][ef_idx] * laser_E0**2)
				#plot
				print(re_plot)
				ax.plot(self.data.smr_lst[0], re_plot,'-', color=dir_color[x]	,label=r'$J_'+dim_str[x]+r'$ @ $E_f=${:+4.2f}'.format(ef_val))
				x_min	=	min(self.data.smr_lst[0])
				x_max	=	max(self.data.smr_lst[0])
				plt.xlabel(r'$\eta$ (eV)', 	fontsize=label_size)
				#~~~
			elif x_axis_para_id==1:		#	VS	FERMI LEVEL
				#collect
				for ef_idx, ef_val in enumerate(self.data.ef_lst[0]):
					re_plot.append(0)
					for i in range(0,dim):
						for j in range(0,dim):
							phi_laser	=	laser.get_phase(i,j)
							re_plot[-1]		= 		re_plot[-1] + np.real(	scale	 *	phi_laser*	self.scndPhoto_data[x][i][j][hw_idx][smr_idx][ef_idx]* laser_E0**2)
				#plot
				ax.plot(self.data.ef_lst[0], re_plot,'-', color=dir_color[x]	,label=r'$J_'+dim_str[x]+r'$ @ $\eta=${:+4.2f}')
				x_min	=	min(self.data.ef_lst[0])
				x_max	=	max(self.data.ef_lst[0])
				plt.xlabel(r'$E_f$ (eV)', 	fontsize=label_size)
				#~~~
			else:						#	VS	FREQUENCY
				#collect
				for hw_idx, hw_val in enumerate(self.data.hw_lst[0]):
					re_plot.append(0)
					for i in range(0,dim):
						for j in range(0,dim):
							phi_laser	=	laser.get_phase(i,j)
							#print(phi_laser)
							re_plot[-1]		= 		re_plot[-1] + np.real(	scale	 *	phi_laser*	self.scndPhoto_data[x][i][j][hw_idx][smr_idx][ef_idx]* laser_E0**2)
				#	plot
				ax.plot(self.data.hw_lst[0], re_plot,'-', color=dir_color[x]	,label=r'$J_'+dim_str[x]+r'$ @ $E_f=${:+4.2f}'.format(ef_val)+r'; $\eta=${:+4.2f}'.format(smr_val))
				#
				x_max	=	max(self.data.hw_lst[0])
				x_min	=	min(self.data.hw_lst[0])
				plt.xlabel(r'$ \hbar \omega $ (eV)', 	fontsize=label_size)
		#~~~~~~~
		#~~~~~~~
		#~~~~~~~
		#~~~~~~~
		#
		# handle x-axis
		ax.set_xlim([x_min, x_max])
		#
		# handle y-axis
		ax.set_ylim([lower_bound, upper_bound  ])
		ax.set(ylabel=r'$J^a{'+'}\;$' +	self.unit_str)
		print("WARNING units are fucked at the moment TODO")
		ax.yaxis.label.set_size(label_size)
		#
		# handle legend
		if plot_legend:
			plt.legend(loc='lower right')

		#
		#handle title
		title = title +r'; $\lambda$='+'{:+1d}'.format(laser_lmbda)
	
		if x_axis_para_id==0 or x_axis_para_id==1:
			title	=	title + r' @ $\hbar \omega=$'+r'{:4.2f}'.format(hw_val)
		plt.title(title)
		#
		# WRITE TO FILE
		outFile_path	= self.plot_dir+'/Jphoto_vs_'+para_str+'.pdf'
		plt.tight_layout()
		plt.savefig(outFile_path)
		plt.close()
		#
		print("-------------------------------------------------------------------------------")
		print("")
		print("para info:")
		if not x_axis_para_id==0:
			print("[plot_photoC]: selected smearing: "	,	self.data.smr_lst[0][smr_idx],		" eV")
		if not x_axis_para_id==1:
			print("[plot_photoC]: fermi_level:"			,	self.data.ef_lst[0][ef_idx],		" eV")
		if x_axis_para_id==0 or x_axis_para_id==1:
			print("[plot_photoC]: selected frequency: "	,	self.data.hw_lst[0][hw_idx],		" eV")
		print("")
		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#	end of class 	HW_probe
#
#**************************************************************************************************************************************
#**************************************************************************************************************************************
#--------------------------------------------------------------------------------------------------------------------------------------
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv



#
#
#^^^^
def plot_scnd_photo(root_dir='.',SI=True):
	#	use the above class in here to plot data in folder root_dir
	print('[plot_scnd_photo]:	hello there')
	if os.path.isdir(root_dir):
		#	read data
		newPlot	= merry_plotter(root_dir,SI=SI)
		#
		frank_I		=	10. 						# G W / cm**2		=	1e9	W/cm**2	= 1e9 1e-4 W/m**2	=	1e5 W/m
		frank_I_SI	=	frank_I * 1e4				# G	W / m**2
		frank_I_SI	=	frank_I_SI * 1e9			#	W / m**2
		frank_LASER	=	laser(x=frank_I_SI,x_is_intensity=True)
		#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		#
		#	PLOT RESPONSES 
		#
		#	~~~~~~~~~~~~~~~~~~~~~~~~
		#
		newPlot.plot_photoC(		scale			=	1		, 
									scale_str		=	r'',
									dim				=	3,
									x_axis_para_id	=	123,
									smr_idx=2, ef_idx=1, hw_idx=15,
									laser=frank_LASER			, 
									laser_dir=2					,
									laser_lmbda=-1,
									line_width=1.5,label_size=14, xtick_size=12, ytick_size=12, marker_size=1.1,
									upper_bound		=	1e12,
									lower_bound		=	-1e12,
									plot_legend=True,
									title='LASER pol x'
							)
		print("...")
		print('[plot_scnd_photo]:	all done by')
		#	~~~~~~~~~~~~~~~~~~~~~~~~
		#
		#----------------------------------------------------------------------------------------------------------------------
	else:
		print('[plot_scnd_photo]:	ERROR '+str(root_dir)+'	seems to be non existing. please specify valid folder')
#
#
plot_scnd_photo(root_dir='.',SI=True)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~










