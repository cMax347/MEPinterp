import os
import sys
import shutil
import math
import numpy as np
from fortran_io		import read_real_tens_file, read_cmplx_tens_file
import matplotlib.pyplot as plt

#***************************************************************************************************************************************************
#***************************************************************************************************************************************************

class conv_data:

	def __init__(self, root_dir):
		self.root_dir		=	root_dir
		self.plot_dir		=	root_dir+'/plots'
		self.dim_string 	= 	['x','y','z']

		if os.path.isdir(self.plot_dir):
			shutil.rmtree(self.plot_dir)
		os.mkdir(self.plot_dir)

		self.dir_lst		= 	[]
		self.nK_lst			=	[]
		#
		self.mep_2014_lst	=	[]
		self.mep_tot_lst 	=	[]
		self.mep_cs_lst	 	=	[]
		self.mep_ic_lst	 	=	[]
		self.mep_lc_lst	 	=	[]
		#
		self.ahc_lst		=	[]
		self.ahc_velo_lst	=	[]
		self.ohc_lst		=	[]
		#
		self.optS_lst		=	[]
		self.optA_lst		=	[]


	
	def collect_data(self):
		for subdir, dirs, files in os.walk(self.root_dir):
			if "nK" in subdir and "w90files" not in subdir and "out" not in subdir and "raw" not in subdir: 
				nK_new	= int(subdir.split("nK")[1])
				print("interpreted nK=",nK_new)
				#
				self.dir_lst.append(subdir)
				self.nK_lst.append(nK_new)
				#
				#self.mep_2014_lst.append(	0.0)
				#self.mep_2014_lst.append(	read_real_tens_file(	subdir + 	'/out/mep/mep_14.dat'			,		'mep'		)			)
				self.mep_tot_lst.append(	read_real_tens_file(	subdir + 	'/out/mep/mep_tens.dat'			, 		'mep'		)			)
				self.mep_cs_lst.append(		read_real_tens_file(	subdir + 	'/out/mep/mep_cs.dat'			, 		'mep'		)			)
				self.mep_ic_lst.append(		read_real_tens_file(	subdir + 	'/out/mep/mep_ic.dat'			, 		'mep'		)			)
				self.mep_lc_lst.append(		read_real_tens_file(	subdir + 	'/out/mep/mep_lc.dat'			, 		'mep'		)			)
				#
				self.ahc_lst.append(		read_real_tens_file(	subdir +	'/out/ahc/ahc_tens.dat'			,		'ahc'		)			)
				self.ahc_velo_lst.append(	read_cmplx_tens_file(	subdir +	'/out/ahc/ahc_velo.dat'			,		'ahcVELO'	)			)
				self.ohc_lst.append(		read_cmplx_tens_file(	subdir + 	'/out/ahc/ohc_kubo.dat'			,		'ohcVELO'	)			)
				#
				self.optS_lst.append(		read_cmplx_tens_file(	subdir + 	'/out/opt/opt_Ssymm.dat'		,		'optS'		)			)
				self.optA_lst.append(		read_cmplx_tens_file(	subdir + 	'/out/opt/opt_Asymm.dat'		,		'optA'		)			)


		#sort by number of kpts used\
		print("raw nK_lst=" + str(self.nK_lst))
		zipped	= zip(self.nK_lst, self.mep_tot_lst, self.mep_cs_lst, self.mep_ic_lst, self.mep_lc_lst, self.ahc_lst)
		

		sort 	= sorted(zipped)
		self.nK_lst,  self.mep_tot_lst, self.mep_cs_lst, self.mep_ic_lst, self.mep_lc_lst, self.ahc_lst	= map(list,zip(*sort))

		print("sorted nK_lst=" + str(self.nK_lst))
		nK_per_dim_lst	=	self.get_nK_plot(False)
		print("sorted nK_per_dim_lst=" + str(nK_per_dim_lst))

		


#***************************************************************************************************************************************************
#***************************************************************************************************************************************************
#			PLOTTING	
#***************************************************************************************************************************************************
#***************************************************************************************************************************************************

	def get_nK_plot(self, show_tot_nK):
		nK_plot	= 	[]
		for nk in self.nK_lst:
			if show_tot_nK:
				nK_plot.append(nk)
			else:
				nK_plot.append(math.ceil(np.power(float(nk),1./3.)))
		#
		return nK_plot




	def plot_mep(self, tick_label_size=12, show_tot_nK=True, show_indi_mep=True):
		#prepare uniform k_plot list
		nK_plot	= 	self.get_nK_plot(show_tot_nK)
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
					mep_14_ab	= []
					mep_py_ab	= []
					for mep_cs in self.mep_cs_lst:
						mep_cs_ab.append(	mep_cs[a][b])	
					for mep_ic in self.mep_ic_lst:
						mep_ic_ab.append(	mep_ic[a][b])	
					for mep_lc in self.mep_lc_lst:
						mep_lc_ab.append(	mep_lc[a][b])	

					#for idx, mep_tot in enumerate(self.mep_cs_lst):
					#	mep_py_ab.append(	self.mep_lc_lst[idx][a][b] + self.mep_ic_lst[idx][a][b] + self.mep_cs_lst[idx][a][b]	)

					#for mep_14 in self.mep_2014_lst:
					#	mep_14_ab.append(	mep_14[a][b])
				#plot
				fig, ax  = plt.subplots(1,1) 
				plt.semilogx(nK_plot, mep_tot_ab, '+-'	,color='black' ,label="tot")
				if show_indi_mep:
					plt.semilogx(nK_plot, mep_cs_ab,	'--',	color='red',		label="CS"	)					
					plt.semilogx(nK_plot, mep_lc_ab,	'o-',	color='darkblue',		label="LC"	)
					plt.semilogx(nK_plot, mep_ic_ab,	'^-',	color='lightblue',		label="IC"	)
					
					#plt.semilogx(nK_plot, mep_py_ab , '+-', color='orange', label="py sumed")
					

					#plt.semilogx(nK_plot, mep_14_ab,	'v-',	color='magenta',	label="2014")
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
				plt.close()
	#***************************************************************************************************************************************************





	def get_min_max_entry(self,	tens_lst	):
		raw_entry_lst	=	[]
		#
		for tens in tens_lst:
			for row in tens:
				for elem in row:
					raw_entry_lst.append(	np.real(elem)	)
					if isinstance(elem,complex):
						raw_entry_lst.append(	np.imag(elem)	)
		#
		min_val	=	min(	raw_entry_lst	)
		max_val	=	max(	raw_entry_lst	)
		#
		return min_val, max_val




	def plot_ahc(self, tick_label_size=12,  show_tot_nK=True, ylim=["null","null"]):
		#prepare uniform k_plot list
		nK_plot	= 	self.get_nK_plot(show_tot_nK)
		##
		#
		min_ahc, max_ahc		=	self.get_min_max_entry(	self.ahc_lst	)
		#
		min_ahcK, max_ahcK		=	self.get_min_max_entry(	self.ahc_velo_lst	)
		#
		min_ohc, max_ohc		=	self.get_min_max_entry(	self.ohc_lst	)
		#
		#
		max_plot	=	max(	[max_ahc, max_ahcK, max_ohc])
		min_plot	=	min(	[min_ahc, min_ahcK, min_ohc])
		#
		##
		print("ahc tens range= ["+str(min_ahc)+", "+str(max_ahc)+"].")
		print("ahc_kubo (optical hall cond. WX) plot range= ["+str(min_ahcK)+", "+str(max_ahcK)+"].")
		print("ohc_kubo (optical conductivity w90) plot range= ["+str(min_ohc)+", "+str(max_ohc)+"].")
		print("ahc plot range= ["+str(min_plot)+", "+str(max_plot)+"].")

		#
		for a in range(0,3):
			for b in range(0,3):
				#
				ahc_ab				= 	[]
				re_ahc_velo_ab		=	[]
				im_ahc_velo_ab		=	[]				
				re_ohc_ab			=	[]
				im_ohc_ab			=	[]

				for ahc_tens in self.ahc_lst:
					ahc_ab.append(		ahc_tens[a][b]		)	
				for ahc_velo in self.ahc_velo_lst:
					re_ahc_velo_ab.append(		np.real(	ahc_velo[a][b]	)		)
					im_ahc_velo_ab.append(		np.imag(	ahc_velo[a][b]	)		)
				for ohc_tens in self.ohc_lst:
					re_ohc_ab.append(			np.real(	ohc_tens[a][b]	)		)
					im_ohc_ab.append(			np.imag(	ohc_tens[a][b]	)		)					
				

				#plot
				fig, ax  = plt.subplots(1,1) 
				plt.semilogx(nK_plot, ahc_ab,	'+-',	color='black', label="AHC" )
				#
				plt.semilogx(nK_plot, re_ahc_velo_ab		,	'^-'	,		color='blue'	,	label='RE OHC')
				plt.semilogx(nK_plot, im_ahc_velo_ab		,	'v-'	,		color='orange'	,	label='IM OHC')
				#
				plt.semilogx(nK_plot, re_ohc_ab				,	'x-'	,		color='red',	label='RE OC')
				plt.semilogx(nK_plot, im_ohc_ab				,	'*-'	,		color='green',	label='IM OC')



				plt.ylim(	[	-.3,	.3	]	)

				#aesthetics
				#plt.xticks(np.array([nK_plot]))

				plt.tick_params(axis='both', which='both',left=True,right=True, direction='in',labelsize=tick_label_size)
				plt.ylabel(r'$\rho$_'+self.dim_string[a]+self.dim_string[b]+' (arb.unit TODO)')
				if show_tot_nK:
					plt.xlabel("nK")
				else:
					plt.xlabel("nK/dim")
				plt.legend()
				plt.title("Hall like Kubo response tensor")
				plt.tight_layout()
				plt.savefig(self.plot_dir+'/ahc_'+self.dim_string[a]+self.dim_string[b]+'_k_conv.pdf')
				#
				print("saved AHC plot")
				plt.close()

	#***************************************************************************************************************************************************


	def plot_opt(self, tick_label_size=12, show_tot_nK=True):
		#prepare uniform k_plot list
		nK_plot	= 	self.get_nK_plot(show_tot_nK)
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
				plt.close()



#***************************************************************************************************************************************************
#***************************************************************************************************************************************************
#***************************************************************************************************************************************************
#***************************************************************************************************************************************************
#***************************************************************************************************************************************************
#***************************************************************************************************************************************************


#default data folder
data_folder		=	"k_conv_cluster_phi0.0"


if len(sys.argv) >	1:
	data_folder	=	sys.argv[1]


if os.path.isdir(data_folder):
	#plot data within folder
	cluster_data	=	conv_data(data_folder)
	cluster_data.collect_data()

	cluster_data.plot_mep( tick_label_size=12,	show_tot_nK=False	, show_indi_mep=True	)
	cluster_data.plot_ahc( tick_label_size=12,	show_tot_nK=False							)
	cluster_data.plot_opt( tick_label_size=12,	show_tot_nK=False							)
else:
	print("ERROR ./"+data_folder+"	is not a directory. Please provide valid folder as cli argument")


