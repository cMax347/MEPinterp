import numpy as np
import datetime
import os
from mep_worker import MEP_worker 
import matplotlib.pyplot as plt





class Phi_probe:

	def __init__(		self,n_phi, val_bands, mp_grid, gamma_scale, kubo_tol, 	
						hw, eFermi, Tkelvin, eta_smearing,debug_mode, do_gauge_trafo='T' ,
						do_write_velo='F', do_write_mep_bands='F',
						do_mep='T', do_kubo='F', do_ahc='F', do_opt='F', do_gyro='F'
						):
		self.n_phi 			= 	n_phi
		self.val_bands		= 	val_bands
		self.mp_grid		= 	mp_grid
		self.gamma_scale	=	gamma_scale
		self.kubo_tol		= 	kubo_tol
		self.hw				= 	hw 
		self.eFermi			= 	eFermi 
		self.Tkelvin		= 	Tkelvin 
		self.eta_smearing	= 	eta_smearing
		self.debug_mode		= 	debug_mode 
		self.do_gauge_trafo	= 	do_gauge_trafo
		self.do_write_velo	=	do_write_velo
		self.do_write_mep_bands	=	do_write_mep_bands
		self.do_mep			=	do_mep
		self.do_kubo		=	do_kubo
		self.do_ahc			=	do_ahc
		self.do_opt			=	do_opt
		self.do_gyro		=	do_gyro

		#derived attributes
		self.root_dir	= os.getcwd()+'/'+datetime.date.today().strftime("%d%B%Y")+'_mp'+str(self.mp_grid[0])+str(self.mp_grid[1])+str(self.mp_grid[2])
		
		#search new root folder name
		if os.path.isdir(self.root_dir):
			old_path 	= self.root_dir
			cnt 		= 0
			while os.path.isdir(old_path) and cnt < 5:
				old_path	= old_path + '.old'
				cnt			= cnt + 1
				try:
					os.rename(self.root_dir, old_path)
				except OSError:
					print(old_path+ ' exists already')

		self.plot_dir		= 	self.root_dir+'/mep_tot_plots'
		self.phi_tot_data 	= 	[]
		self.phi_cs_data 	= 	[]
		self.phi_lc_data 	= 	[]
		self.phi_ic_data 	= 	[]
		self.phi_bands_data	=	[]




	def __del__(self):
		print('done with Phi probing, by\n\n')





	def iterate_phi(self,mpi_np=1, plot_bandstruct=False):
		for phi  in np.linspace(0.0, 2.0, num = self.n_phi):		#iterate over relative phi (phi_rel = phi / np.pi)
			work_dir =	self.root_dir+'/phi'+str(phi)		
			phi_pi	 = 	phi
			#init current job
			worker = MEP_worker(	self.root_dir, 
									work_dir, 
									phi_pi, 						#give  phi_rel*pi to calculation
									self.val_bands, 
									self.mp_grid, 
									self.kubo_tol,
									self.hw,
									self.eFermi,
									self.Tkelvin,
									self.eta_smearing,
									self.debug_mode,
									self.do_gauge_trafo,
									self.do_write_velo,
									self.do_write_mep_bands,
									self.do_mep,
									self.do_kubo,
									self.do_ahc,
									self.do_opt,
									self.do_gyro
								)
			#run calc 
			worker.run(mpi_np=mpi_np)
			mep_tens, mep_cs, mep_lc, mep_ic, mep_bands = worker.get_mep_tens()
			self.phi_tot_data.append(		[phi, mep_tens	]			)
			self.phi_cs_data.append(		[phi, mep_cs	]			)
			self.phi_lc_data.append(		[phi, mep_lc	]			)
			self.phi_ic_data.append(		[phi, mep_ic	]			)
			self.phi_bands_data.append(		[phi, mep_bands ]			)
			#plot bands
			if plot_bandstruct:
				try:
					worker.plot_bands()
				except:
					print('could not plot bands')
			
		self.phi_tot_data 	= 	sorted(self.phi_tot_data)
		self.phi_cs_data	=	sorted(self.phi_cs_data)
		self.phi_lc_data	=	sorted(self.phi_lc_data)
		self.phi_ic_data	=	sorted(self.phi_ic_data)
		self.phi_bands_data	=	sorted(self.phi_bands_data)



	def print_results_container(self):
		print('^^^^^raw data ^^^^^^^')
		print('mep_bands :')
		print(self.phi_bands_data)
		#print(self.phi_bands_data[0][0])

		print('phi_tot_data =',self.phi_tot_data)
		print('------------------------------------')




	

	

	def plot_mep_over_phi(self, plot_contributions=True, plot_band_res=False, label_size=14, xtick_size=12, ytick_size=12):
		try:
			os.mkdir(self.plot_dir)
		except OSError:
			print('Could not make directory ',self.plot_dir)

		phi_plot 		= []
		phi_plot_bands	= []
		mep_tot_data 	= []
		mep_cs_data		= []
		mep_lc_data		= []
		mep_ic_data		= []
		mep_band_data	= []
		for data in self.phi_tot_data:
			phi_plot.append(data[0])
			mep_tot_data.append(data[1])
	
		if plot_band_res:
			for n_phi,mep_band in enumerate(self.phi_bands_data):
				print('current PHI: #',n_phi)
				phi_plot_bands.append(mep_band[0])
				mep_band_data.append(mep_band[1])


		print('phi_plot=',phi_plot)
		print('mep_tot_data=',mep_tot_data)
		
		if plot_band_res:
			print('phi_plot_bands=',phi_plot_bands)
			print('mep_band_data=',mep_band_data)


		if plot_contributions:
			for data in self.phi_cs_data:
				mep_cs_data.append(data[1])
			for data in self.phi_lc_data:
				mep_lc_data.append(data[1])
			for data in self.phi_ic_data:
				mep_ic_data.append(data[1])


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
				mep_band_plot	= []
				mep_cs_plot		= []
				mep_lc_plot		= []
				mep_ic_plot		= []
				for mep_tens in mep_tot_data:
					mep_tot_plot.append(				self.gamma_scale *		mep_tens[i][j]	)
					

				#collect band resolved data	
				if plot_band_res:
					print('i=',i,' j=',j, '	band contributions:')
					for phi, data in enumerate(mep_band_data):
						#print('mep_band_data: #phi=',phi)
						mep_band_plot.append([])
						print('---')
						n_bands	= 0
						for idx,mep_band in enumerate(data):
							n_bands	=	n_bands	+ 1
							#print('band #',idx,'	data:',mep_band)
							mep_band_plot[-1].append(	mep_band[i][j])
							print('phi=',phi,' #band',idx, ' mep_band: ',mep_band[i][j])
					print('------')

					print("#bands	=	"+str(n_bands))
					print("mep_band_plot:",mep_band_plot)
					mep_band_final	= []

					#n_bands	= 2
					for band in range(n_bands):
						mep_band_final.append([])
						for phi in phi_plot_bands:
							mep_band_final[-1].append(phi)

					print('mep_band_final container:',	mep_band_final)
					
					for phi,mep_bands in enumerate(	mep_band_plot	):
						print("phi=",phi,	'mep:')
						print(mep_bands)
						for band, mep in enumerate(	mep_bands):
							mep_band_final[band][phi]	=	mep

					print('mep_band_final values:',	mep_band_final)

				if plot_contributions:
					for mep_cs in mep_cs_data:
						mep_cs_plot.append(				self.gamma_scale *		mep_cs[i][j]	)
					for mep_lc in mep_lc_data:
						mep_lc_plot.append(				self.gamma_scale *		mep_lc[i][j]	)
					for mep_ic in mep_ic_data:
						mep_ic_plot.append(				self.gamma_scale *		mep_ic[i][j]	)


				#print('mep_tot_plot=',mep_tot_plot)
				#print('mep_band_plot=',mep_band_plot)



				#do PLOT
				fig, ax  = plt.subplots(1,1) 
				plt.plot(phi_plot, mep_tot_plot,'-+', color='black',label="tot")

				if plot_contributions:
					plt.plot(phi_plot, mep_cs_plot,		'--', 	color='red',		label="CS")
					plt.plot(phi_plot, mep_lc_plot,		'o-', 	color='darkblue',	label="LC")
					plt.plot(phi_plot, mep_ic_plot,		'^-', 	color='lightblue',	label="IC")

				if plot_band_res:
					for band in range(n_bands):
						plt.plot(phi_plot_bands,	mep_band_final[band],	'x-', label=" #band "+str(band)		)

				
				#X-AXIS
				plt.xlabel(r'$\varphi$',	fontsize=label_size)
				ax.set_xlim([0,2])
				ax.set_xticks(np.array([0,0.5,1,1.5,2]))
				ax.set_xticklabels(np.array([r'$0$','',r'$\pi$','',r'$2\pi$']))
				plt.tick_params(axis='x',which='major', direction='in',labelsize=xtick_size)
				
				#Y-AXIS
				plt.ylabel(r'$\alpha_{'+dim_str[i]+dim_str[j]+'}$',	fontsize=label_size)
				ax.set_ylim([self.gamma_scale *  mep_min,self.gamma_scale *  mep_max])
				plt.tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)

				if plot_contributions or plot_band_res:
					plt.legend()

				plt.title('gamma_scale='+str(self.gamma_scale))

				plt.tight_layout()
				plt.savefig(self.plot_dir+'/mep_'+dim_str[i]+dim_str[j]+'.pdf')
				plt.close()


				

	








def probe_phi(		n_phi, val_bands, mp_grid,mpi_np=1, gamma_scale=1,
					kubo_tol=1e-3, hw=0.001, eFermi=0.0, Tkelvin=11.0, eta_smearing=0.2, 
					plot_bandstruct	=	False, 
					plot_orb_cont	=	True,
					plot_band_res	=	False,
					debug_mode		=	True, 
					do_gauge_trafo	=	True,
					do_write_velo	=	False,
					do_write_mep_bands	= False,
					do_mep			=	'T', 
					do_kubo			=	'F', 
					do_ahc			=	'F', 
					do_opt			=	'F', 
					do_gyro			=	'F'	
				):
	myTest	= Phi_probe(	n_phi, val_bands, mp_grid, gamma_scale, 
							kubo_tol, hw, eFermi, Tkelvin, eta_smearing,
							debug_mode, do_gauge_trafo, do_write_velo, do_write_mep_bands,
							do_mep, do_kubo, do_ahc, do_opt, do_gyro	
						)
	#
	myTest.iterate_phi(plot_bandstruct=plot_bandstruct, mpi_np=mpi_np)
	myTest.print_results_container()
	#try:
	#myTest.plot_mep_over_phi(	plot_contributions	=	plot_orb_cont,	
	#							plot_band_res		= 	plot_band_res, 
	#							label_size			=	14, 
	#							xtick_size			=	12, 
	#							ytick_size			=	12
	#						)
	#except:
	#	print("plotting failed. Please try plotting with plot_probe_phi.py")
	#finally:
	#	print('')
	#	print('')
	#	print('all done')






probe_phi(		#parameter space probe density:
				n_phi				=	21					, 
				val_bands			=	2					, 
				#numerical parameters
				mp_grid				=	[16,16,16]			, 
				mpi_np				=	4					,
				gamma_scale			=	1.0					,
				kubo_tol			=	1e-5				, 
				hw					=	0.0					, 
				eFermi				=	0.0					, 
				Tkelvin				=	10.0				,
				eta_smearing		=	0.1					, 
				#set properties of plots
				plot_bandstruct		=	False				, 
				plot_orb_cont		=	True				,				
				plot_band_res		=	False				,
				#additional fortran controllers
				debug_mode			=	True				, 
				do_gauge_trafo		=	True				,	
				do_write_velo		=	False				,
				do_write_mep_bands	=	True				,
				do_mep				=	'T'					, 
				do_kubo				=	'F'					, 
				do_ahc				=	'F'					, 
				do_opt				=	'F'					, 
				do_gyro				=	'F'
			)



