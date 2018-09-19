import os
import shutil
import numpy as np
from fortran_io		import read_real_tens_file, read_cmplx_tens_file
import matplotlib.pyplot as plt

#***************************************************************************************************************************************************
#***************************************************************************************************************************************************

class conv_data:

	def __init__(self, root_dir):
		self.root_dir	=	root_dir
		self.plot_dir	=	root_dir+'/plots'
		self.dim_string = 	['x','y','z']

		if os.path.isdir(self.plot_dir):
			shutil.rmtree(self.plot_dir)
		os.mkdir(self.plot_dir)

		self.dir_lst	= 	[]
		self.nK_lst		=	[]
		#
		self.mep_tot_lst=	[]
		self.mep_cs_lst	=	[]
		self.mep_ic_lst	=	[]
		self.mep_lc_lst	=	[]
		#
		self.ahc_lst	=	[]
		self.optS_lst	=	[]
		self.optA_lst	=	[]


	
	def collect_data(self):
		for subdir, dirs, files in os.walk(self.root_dir):
			if "nK" in subdir and "w90files" not in subdir and "out" not in subdir and "raw" not in subdir: 
				nK_new	= int(subdir.split("nK")[1])
				print("interpreted nK=",nK_new)
				#
				self.dir_lst.append(subdir)
				self.nK_lst.append(nK_new)
				#
				self.mep_tot_lst.append(	read_real_tens_file(subdir + '/out/mep/mep_tens.dat'		, 		'mep')			)
				self.mep_cs_lst.append(		read_real_tens_file(subdir + '/out/mep/mep_cs.dat'			, 		'mep')			)
				self.mep_ic_lst.append(		read_real_tens_file(subdir + '/out/mep/mep_ic.dat'			, 		'mep')			)
				self.mep_lc_lst.append(		read_real_tens_file(subdir + '/out/mep/mep_lc.dat'			, 		'mep')			)
				#
				self.ahc_lst.append(		read_real_tens_file(subdir +	'/out/ahc/ahc_tens.dat'		,		'ahc')			)
				#
				self.optS_lst.append(		read_cmplx_tens_file(subdir + '/out/opt/opt_Ssymm.dat'		,		'optS')			)
				self.optA_lst.append(		read_cmplx_tens_file(subdir + '/out/opt/opt_Asymm.dat'		,		'optA')			)


		#sort by number of kpts used\
		print("raw nK_lst=" + str(self.nK_lst))
		zipped	= zip(self.nK_lst, self.mep_tot_lst, self.mep_cs_lst, self.mep_ic_lst, self.mep_cs_lst, self.ahc_lst)
		sort 	= sorted(zipped)
		self.nK_lst, self.mep_tot_lst, self.mep_cs_lst, self.mep_ic_lst, self.mep_cs_lst, self.ahc_lst	= map(list,zip(*sort))
		print("sorted nK_lst=" + str(self.nK_lst))

		


#***************************************************************************************************************************************************
#***************************************************************************************************************************************************
#			PLOTTING	
#***************************************************************************************************************************************************
#***************************************************************************************************************************************************

	def plot_mep(self, tick_label_size=12, show_indi_mep=True, show_tot_nK=True):
		#prepare uniform k_plot list
		nK_plot	= 	[]
		for nk in self.nK_lst:
			if show_tot_nK:
				nK_plot.append(nk**3)
			else:
				nK_plot.append(nk)
		#
		for a in range(0,3):
			for b in range(0,3):
				#
				mep_tot_ab	= []
				for mep_tens in self.mep_tot_lst:
					mep_tot_ab.append(	mep_tens[a][b]		)
				if show_indi_mep:
					mep_cs_ab	= []
					mep_ic_ab	= []
					mep_lc_ab	= []
					for mep_cs in self.mep_cs_lst:
						mep_cs_ab.append(	mep_cs[a][b])	
					for mep_ic in self.mep_ic_lst:
						mep_ic_ab.append(	mep_ic[a][b])	
					for mep_lc in self.mep_lc_lst:
						mep_lc_ab.append(	mep_lc[a][b])	
				#plot
				fig, ax  = plt.subplots(1,1) 
				plt.semilogx(nK_plot, mep_tot_ab, '+-'	,color='black' ,label="tot")
				if show_indi_mep:
					plt.semilogx(nK_plot, mep_cs_ab,	'--',	color='red',		label="CS"	)					
					plt.semilogx(nK_plot, mep_lc_ab,	'o-',	color='darkblue',		label="LC"	)
					plt.semilogx(nK_plot, mep_ic_ab,	'^-',	color='lightblue',		label="IC"	)
				#aesthetics
				plt.tick_params(axis='both', which='both',left=True,right=True, direction='in',labelsize=tick_label_size)
				plt.ylabel(r'$\alpha$_'+self.dim_string[a]+self.dim_string[b]+' (a.u.)')
				if show_tot_nK:
					plt.xlabel("nK")
				else:
					plt.xlabel("nK/dim")
				plt.legend()
				plt.title("MEP response tensor")
				plt.tight_layout()
				plt.savefig(self.plot_dir+'/mep_'+self.dim_string[a]+self.dim_string[b]+'_k_conv.pdf')
				#
				print("saved MEP plot")
	#***************************************************************************************************************************************************

	def plot_ahc(self, tick_label_size=12,  show_tot_nK=True):
		#prepare uniform k_plot list
		nK_plot	= 	[]
		for nk in self.nK_lst:
			if show_tot_nK:
				nK_plot.append(nk**3)
			else:
				nK_plot.append(nk)
		#
		for a in range(0,3):
			for b in range(0,3):
				#
				ahc_ab		= []
				for ahc_tens in self.ahc_lst:
					ahc_ab.append(		ahc_tens[a][b]		)	
				#plot
				fig, ax  = plt.subplots(1,1) 
				plt.semilogx(nK_plot, ahc_ab,	'+-',	color='black', label="AHC" )
				#aesthetics
				plt.tick_params(axis='both', which='both',left=True,right=True, direction='in',labelsize=tick_label_size)
				plt.ylabel(r'$\rho$_'+self.dim_string[a]+self.dim_string[b]+' (arb.unit TODO)')
				if show_tot_nK:
					plt.xlabel("nK")
				else:
					plt.xlabel("nK/dim")
				plt.legend()
				plt.title("AHC response tensor")
				plt.tight_layout()
				plt.savefig(self.plot_dir+'/ahc_'+self.dim_string[a]+self.dim_string[b]+'_k_conv.pdf')
				#
				print("saved AHC plot")

	#***************************************************************************************************************************************************


	def plot_opt(self, tick_label_size=12, show_tot_nK=True):
		#prepare uniform k_plot list
		nK_plot	= 	[]
		for nk in self.nK_lst:
			if show_tot_nK:
				nK_plot.append(nk**3)
			else:
				nK_plot.append(nk)
		#
		for a in range(0,3):
			for b in range(0,3):
				#grep component ab
				re_opt_s_ab	= []
				im_opt_s_ab	= []
				re_opt_a_ab	= []
				im_opt_a_ab = []
				#
				for optS_tens in self.optS_lst:
					re_opt_s_ab.append(		np.real(	optS_tens[a][b]	))
					im_opt_s_ab.append(		np.imag(	optS_tens[a][b]	))
				for optA_tens in self.optA_lst:
					re_opt_a_ab.append(		np.real(	optA_tens[a][b]	))
					im_opt_a_ab.append(		np.imag(	optA_tens[a][b]	))

			
				fig, ax	=	plt.subplots(1,1)
				plt.semilogx(nK_plot, re_opt_s_ab, '^-',color="blue", label='Re symm')
				plt.semilogx(nK_plot, im_opt_s_ab, 'v-',color="orange" , label='Im symm')
				plt.semilogx(nK_plot, re_opt_s_ab, 'o-',color="green", label='Re asymm')
				plt.semilogx(nK_plot, re_opt_s_ab, 's-',color="lime", label='Im asymm')
				#aesthetics
				plt.tick_params(axis='both', which='both',left=True,right=True, direction='in',labelsize=tick_label_size)
				plt.ylabel(r'$\rho$_'+self.dim_string[a]+self.dim_string[b]+' (arb.unit TODO)')
				if show_tot_nK:
					plt.xlabel("nK")
				else:
					plt.xlabel("nK/dim")
				plt.legend()
				plt.title("optical conductivity tensor")
				plt.tight_layout()
				plt.savefig(self.plot_dir+'/opt_'+self.dim_string[a]+self.dim_string[b]+'_k_conv.pdf')
				#
				print("saved OPT plot")



#***************************************************************************************************************************************************
#***************************************************************************************************************************************************
#***************************************************************************************************************************************************
#***************************************************************************************************************************************************
#***************************************************************************************************************************************************
#***************************************************************************************************************************************************



cluster_data	=	conv_data("k_conv_cluster_phi0.0")
cluster_data.collect_data()
cluster_data.plot_mep( tick_label_size=12, show_indi_mep=True, show_tot_nK=True)
cluster_data.plot_ahc( tick_label_size=12, show_tot_nK=False)
cluster_data.plot_opt( tick_label_size=12, show_tot_nK=False)


