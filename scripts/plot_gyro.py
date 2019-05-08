import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib2tikz


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#	constants etc.
#____
gyroD_id	=	'gyroD'
dim_str		=	['x','y','z']
#
cond_quantum		=	2.434135 * 1e-4		#	1	[	e**2/hbar	]_atomic	=	2.434135×10^-4 	[	S	]_SI
elem_e_over_hartree	=	0.03674932			#	1	[	e/E_h		]_atomic	=	0.03674932 		[	1/V	]_SI
#
omega_au_to_si		=	cond_quantum	*	elem_e_over_hartree * 1e6 		#[	e**2/hbar e/E_h	]_atomic	-> [1e-6 A/V**2] = [mu A/V**2] 
#
hartree_to_ev		=	27.211385			#	1 [E_h]_atomic	= 27.211385 [eV]

print("[plot_gyro]: WARNING cond_quantum not set correctly ")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#	functions
#____
def my_levi_civ(i,j,k):
	return int((i-j)*(j-k)*(k-i)/2)


def cvrt_Dtens_to_nonlin_AHE(D_tens, tau_au=1. ,omega_au=0):
	#	see eq.(B7) in	Souza et al., PRB 97, 035158 (2018)
	#		
	#		=> omega	=
	omega	=	np.zeros((3,3,3))		*1j
	#
	tau_cmplx	=	tau_au / (1. - 1j * tau_au * omega_au )

	pre_fact	=	0j	-	 tau_cmplx 
	#	iterate directions	
	for c in range(3):
		for b in range(3):
			for a in range(3):
				#
				#print(	'({:2d} {:2d} {:2d}):'.format(a, b, c))
				for d in range(3):
					omega[a][b][c]	=	omega[a][b][c]		+	my_levi_civ(a+1,d+1,c+1)	*	D_tens[b][d]
					#print(	  '		-> '+str(my_levi_civ(a+1,d+1,c+1))	+ '		->	'+str(my_levi_civ(a+1,d+1,c+1)	*	D_tens[b][d])	)
	#
	return pre_fact*omega



def test_levi_civ():
	print("	BEGIN LEVI CIV TEST:")

	for i in 1, 2, 3:
		for j in 1, 2, 3:
			for k in 1, 2, 3:
				print(	'({:2d} {:2d} {:2d})'.format(i, j, k)	+  '-> {:2d}'.format(my_levi_civ(i, j, k))	)
	print("end levi-civ test")

def print_tens(tens):
	for dim in range(3):
		print(dim_str[dim])
		for x in range(3):
			tmp = []
			for elem in tens[dim][x]:
				tmp.append(np.imag(elem))	
			tmp = tens[dim][x]	
			print('		{:4.2e} {:4.2e}\t{:4.2e} {:4.2e}\t{:4.2e} {:4.2e}'.format(	np.real(tmp[0]),	np.imag(tmp[0]),
																					np.real(tmp[1]),	np.imag(tmp[1]),
																					np.real(tmp[2]),	np.imag(tmp[2])
																					)
				)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def safe_mkdir(new_dir):
	if not os.path.isdir(new_dir):
		try:
			os.mkdir(new_dir)
		except OSError:
			print('[plot_hall_like]:	Could not make directory ',new_dir, '	(proably exists already)')
		finally:
			print("~")
	else:
		print('[plot_hall_like]: '+new_dir+"	exists already! (WARNING older plots might be overwriten)")


def get_tau_au(delta_au):
	#	in general:
	#			tau	=	hbar / omega		
	#	in au (hbar=1):
	#			tau	=	1 	/	omega		
	tau_au	=	1./delta_au
	return tau_au


def loop_freq_si(D_tens, eta_smr_ev=0.01,	freq_min=0, freq_max=6, n_freq=121):
	freq_ev_lst= np.linspace(freq_min,freq_max,n_freq)	
	eta_smr_au	=	eta_smr_ev / hartree_to_ev
	
	om_ahe_si	=	[]
	tau_au		=	get_tau_au(eta_smr_ev/ hartree_to_ev)
	print('realaxation time tau=',tau_au,'x [2.418884326505(16)×10−17 s]')
	for freq_ev in freq_ev_lst:
		
		print(	"hw=",freq_ev,' (eV)')
		#
		ahe_hw_au	=	cvrt_Dtens_to_nonlin_AHE(D_tens,	tau_au	=	tau_au	,	omega_au=	freq_ev / hartree_to_ev )
		om_ahe_si.append(	ahe_hw_au * omega_au_to_si	)

	return freq_ev_lst, om_ahe_si



def plot_hall_like(		 hw_lst,	tens_lst	,plot_dir='plots',	 scale=1.0, 
						line_width=1, label_size=14, xtick_size=12, ytick_size=12,
						marker_size=12,
						re_bound=1, im_bound=1
				):
	print("^")
	print("^")
	print("-------------------------------------------------------------------------------")	
	print("		PLOT NONLINEAR AHE (BASED ON GYROTROPIC FORMALISM)")
	print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print("[plot_hall_like]:	try to create folder where plots should go")
	#
	safe_mkdir(plot_dir)
	#
	hw_max	=	max(hw_lst)
	hw_min	=	min(hw_lst)
	#
	#LOOP SPACIAL COMPONENTS OF TENSOR
	for x in range(0,3):
		for i in range(0,3):
			for j in range(0,3):		
				#COLLECT a_ij(phi=0:n_phi)
				#do PLOT
				fig, ax  = plt.subplots(2,1, sharex=True)
				#
				#
				print('ax:',ax)
				ax = ax.flatten()
				print('flat ax:',ax)
				#
				#	TITLE
				title = ' nonlin AHE'
				fig.suptitle(title, fontsize=16)
				#
				#	COLLECT TENSOR(x,i,j) at different frequencies
				real_part	=	[]
				cplx_part	=	[]
				for hw_idx, tens in enumerate(tens_lst):
					real_part.append(np.real(tens[x][i][j]))
					cplx_part.append(np.imag(tens[x][i][j]))				
				#
				#	PLOT 
				ax[0].plot(hw_lst, real_part,	'-', color='blue', label=r'Re $\sigma^'+dim_str[x]+'_'+dim_str[i]+dim_str[j] +'$'	)#,label=r'$\Re$')
				ax[1].plot(hw_lst, cplx_part,	'-', color='red',label=r'Im $\sigma^'+dim_str[x]+'_'+dim_str[i]+dim_str[j]+'$' 	)#,label=r'$\Im$')
				#					
				#	X-AXIS
				ax[1].set_xticks(np.arange(hw_min, hw_max+1.0, 0.1), minor=True)
				ax[1].set_xlim([hw_min, hw_max])
				plt.xlabel(r'$ \hbar \omega $ (eV)',	fontsize=label_size)
				ax[0].set_ylabel(r'$ \sigma^\Re  \; (\mu \mathrm{A}/\mathrm{V}^2)$ ',	fontsize=label_size)
				ax[1].set_ylabel(r'$ \sigma^\Im  \; (\mu \mathrm{A}/\mathrm{V}^2)$ ',	fontsize=label_size)
				#
				#	AXIS GRID?
				#ax.set_ylim([mep_min,mep_max])
				#ax[0].tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)
				#ax[0].tick_params(axis='x',which='both', direction='in',labelsize=xtick_size)
				#ax[1].tick_params(axis='x',which='both', direction='in',labelsize=xtick_size, top=True)
				#ax[1].tick_params(axis='y',which='major', direction='in',labelsize=ytick_size)
				#
				#	LEGENDS
				#handles, labels =	ax.get_legend_handles_labels()
				#fig.legend(handles, labels, loc='lower right')
				#
				#	LAYOUT
				#plt.tight_layout()
				fig.subplots_adjust(hspace=0.01)
				#
				#	WRITE TO FILE
				pdf_outFile_path	= plot_dir+'/Jphoto^'+dim_str[x]+'_'+dim_str[i]+dim_str[j]+'.pdf'
				tikz_outFile_path	= plot_dir+'/Jphoto_'+dim_str[x]+dim_str[i]+dim_str[j]+'.tex'
				plt.savefig(pdf_outFile_path)
				try:
					matplotlib2tikz.save(tikz_outFile_path)
					print("[plot_hall_like]: 	SUCCESS wrote "+tikz_outFile_path)
				except:
					print("[plot_hall_like]: 	WARNING	writing tikz file failed")


				plt.close()
				#
				print('[plot_hall_like]:	finished processing '+dim_str[x]+' -	'+dim_str[i]+dim_str[j]+' tensor, plot saved to: '+pdf_outFile_path	)
	print("-------------------------------------------------------------------------------")
	print("")
	print("")	




#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#	main body
#____
def main():
	#
	#test_levi_civ()
	#
	#	get D_tens
	new_file	=	'./out/gyro_D.npy'
	raw_D_tens	=	np.load(new_file)
	#D_tens		=	read_cmplx_tens_file(new_file,	gyroD_id	)
	#
	print('[plot_gyro]:	raw D_tens:')
	print(raw_D_tens,'\n')
	#
	dimX, dimY, n_ef =	raw_D_tens.shape
	print('input interpretation:	dimX=',dimX)
	print('input interpretation:	dimY=',dimY)
	print('input interpretation:	n_ef=',n_ef)


	D_select	=	raw_D_tens[:,:,0]
	print(D_select.shape)
	#	atomic
	#print('[plot_gyro]:	final nonlinear AHE conductivity (atomic units):')	

	hw_ev_lst, ahe_hw_lst_si	=	loop_freq_si(D_select, eta_smr_ev=0.04,freq_min=0, freq_max=0.6, n_freq=700 )
	
	#
	#	SI 
	print('[plot_gyro]:	final nonlinear AHE conductivity (SI units: mu A/V**2 ):')	
	#
	#	print results to cli 
	#for hw_idx,ahe_si in enumerate(ahe_hw_lst_si):
	#	print('hw=',hw_ev_lst[hw_idx])
	#	print_tens(ahe_si)

	#
	# 	plot responses
	plot_hall_like(	hw_ev_lst,	ahe_hw_lst_si,	 
					plot_dir='./plots', scale=1.0, 
					line_width=1, label_size=14, xtick_size=12, ytick_size=12,
					marker_size=12,
					re_bound=1, im_bound=1
				)


main()