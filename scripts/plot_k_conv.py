import os
import shutil
from fortran_io		import read_mep_file
import matplotlib.pyplot as plt


class conv_data:

	def __init__(self, root_dir):
		self.root_dir	=	root_dir
		self.plot_dir	=	root_dir+'/plots'

		if os.path.isdir(self.plot_dir):
			shutil.rmtree(self.plot_dir)
		os.mkdir(self.plot_dir)

		self.dir_lst	= 	[]
		self.nK_lst		=	[]
		self.mep_tot_lst=	[]
		self.mep_cs_lst	=	[]
		self.mep_ic_lst	=	[]
		self.mep_lc_lst	=	[]


	
	def collect_mep(self):

		for subdir, dirs, files in os.walk(self.root_dir):
			if "nK" in subdir and "w90files" not in subdir and "out" not in subdir and "raw" not in subdir: 
				nK_new	= int(subdir.split("nK")[1])
				print("interpreted nK=",nK_new)

				self.dir_lst.append(subdir)
				self.nK_lst.append(nK_new)

				self.mep_tot_lst.append(	read_mep_file(subdir + '/out/mep_'+		'tens'		+'.dat')			)
				self.mep_cs_lst.append(		read_mep_file(subdir + '/out/mep_'+		 'cs'		+'.dat')			)
				self.mep_ic_lst.append(		read_mep_file(subdir + '/out/mep_'+		 'ic'		+'.dat')			)
				self.mep_lc_lst.append(		read_mep_file(subdir + '/out/mep_'+		 'lc'		+'.dat')			)

		#sort by number of kpts used\
		print("raw nK_lst=" + str(self.nK_lst))
		zipped	= zip(self.nK_lst, self.mep_tot_lst, self.mep_cs_lst, self.mep_ic_lst, self.mep_cs_lst)
		sort 	= sorted(zipped)
		self.nK_lst, self.mep_tot_lst, self.mep_cs_lst, self.mep_ic_lst, self.mep_cs_lst	= map(list,zip(*sort))
		print("sorted nK_lst=" + str(self.nK_lst))

	def plot_conv(self, tick_label_size=12, show_indi_mep=True, show_tot_nK=True):
		dim_string = ['x','y','z']

		#get kpts
		nK_plot	=	[]
		for nk in self.nK_lst:
			if show_tot_nK:
				nK_plot.append(nk**3)
			else:
				nK_plot.append(nk)


		for a in range(0,3):
			for b in range(0,3):
				#grep component ab
				mep_tot_ab	= []
				for mep_tens in self.mep_tot_lst:
					mep_tot_ab.append(	mep_tens[a][b])

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
				plt.ylabel(r'$\alpha$_'+dim_string[a]+dim_string[b]+' (a.u.)')
				if show_tot_nK:
					plt.xlabel("nK")
				else:
					plt.xlabel("nK/dim")
				plt.legend()
				plt.title("MEP response tensor")
				plt.tight_layout()
				plt.savefig(self.plot_dir+'/mep_'+dim_string[a]+dim_string[b]+'_k_conv.pdf')






cluster_data	=	conv_data("k_conv_cluster_phi0.0")
cluster_data.collect_mep()
cluster_data.plot_conv( tick_label_size=12, show_indi_mep=True, show_tot_nK=True)


