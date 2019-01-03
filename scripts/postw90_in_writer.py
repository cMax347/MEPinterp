import datetime


class postw90_job:

	def __init__(self, w90_dir, seedname, val_bands, mp_grid,  hw=0.0,eFermi=0.0, Tkelvin=0.0,eta_smearing=0.0	):
		self.w90_dir	=	w90_dir
		self.seedname	=	seedname
		self.val_bands	=	val_bands
		self.mp_grid	=	mp_grid
		self.hw			=	hw
		self.eFermi		=	eFermi
		self.Tkelvin	=	Tkelvin
		self.eta_smearing=	eta_smearing


	def write_win_file(self ):
		with open(self.w90_dir+'/'+self.seedname+'.win','w') as outfile:
			outfile.write('# input file for postw90 calculation '+'generated on '+datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")+'\n')
			outfile.write('\n')
			#
			outfile.write('berry = true \n'																)
			outfile.write('berry_task = ahc, kubo'														)
			outfile.write('berry_kmesh = '+str(self.mp_grid[0])+' '+str(self.mp_grid[1])+' '+str(self.mp_grid[2])		)

			outfile.write('fermi_energy = '+str(self.eFermi))


			#todo:		-	temperature	(is this via min max value of fermi energy?)
			#			-	fermi smearing for optical conductivity
			#		look how these smearings are implemented in w90



			#
			outfile.write('kubo_freq_max = '+str(self.hw))
			

	









