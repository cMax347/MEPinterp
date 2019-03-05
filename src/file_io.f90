module file_io
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
	use matrix_math,					only:		is_equal_vect, is_herm_mat	
	use constants,						only:		dp, fp_acc, aUtoAngstrm, aUtoEv
	use mpi_community,					only:		mpi_id, mpi_root_id, mpi_nProcs
	use input_paras,					only:		w90_dir, 							&
													out_dir,							& 
													eig_out_dir,						&
													velo_out_dir,						&
													mep_out_dir,						&
													ahc_out_dir,						&
													opt_out_dir,						&
													gyro_out_dir,						&
													raw_dir, 							&
													a_latt, 							&
													debug_mode,							&
													use_R_float,						&
													plot_bands,  phi_laser										

	implicit none

	private
	public									::		read_tb_basis,						& 
													read_kptsgen_pl_file,				&
													!-----------------------------------!
													write_ham_binary,					&
													write_eig_binary,					&
													write_en_binary, read_en_binary,	&
													write_en_global,					&
													write_velo,							&
													!-----------------------------------!
													write_hw_list,						&
													write_mep_bands,					&
													write_mep_tensors,					&
													write_kubo_mep_tensors,				&
													write_ahc_tensor,					&
													write_opt_tensors,					&
													write_gyro_tensors




	interface write_tens_file
		module procedure d2_write_tens_file
		module procedure d3_write_tens_file
		!
		module procedure z_write_tens_file
	end interface write_tens_file



	character(len=64)						::		format='(a,i7.7)'
	integer									::		num_bands
	contains





!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!			WRITE DATA
!--------------------------------------------------------------------------------------------------------------------------------		
!
!
	subroutine write_en_binary(kpt_idx, e_bands)
		integer,		intent(in)		::	kpt_idx
		real(dp),		intent(in)		::	e_bands(:)
		character(len=24)				::	filename
		integer							::	mpi_unit
		!
		mpi_unit	=	100 + mpi_id + 6 * mpi_nProcs
		!
		write(filename, format) raw_dir//'enK.',kpt_idx
		open(unit=mpi_unit,	file = filename, form='unformatted', action='write', access='stream',	status='replace'		)
		write(mpi_unit)	e_bands(:)
		close(mpi_unit) 
		!
		return
	end subroutine


	subroutine write_eig_binary(kpt_idx, U_k)
		integer,		intent(in)		::	kpt_idx
		complex(dp),	intent(in)		::	U_k(:,:)
		character(len=24)				::	filename
		integer							::	mpi_unit, n, m
		!
		mpi_unit	=	100 + mpi_id + 15 * mpi_nProcs
		write(filename, format) eig_out_dir//'eig.',kpt_idx
		!
		open(unit=mpi_unit, file=filename, form='formatted', action='write', access='stream', status='replace')
			do m =1, size(U_k,2)
				do n =1, size(U_k,1)
					write(mpi_unit,'(i2,a,i2,a,e16.7,a,e16.7)') n," ",m," ",real(U_k(n,m),dp)," ",aimag(U_k(n,m))
				end do
			end do
		close(mpi_unit)
		!
		return
	end subroutine


	subroutine write_ham_binary(kpt_idx, H_k)
		integer,		intent(in)		::	kpt_idx
		complex(dp),	intent(in)		::	H_k(:,:)
		character(len=24)				::	filename
		integer							::	mpi_unit, n, m
		!
		mpi_unit	=	100 + mpi_id + 15 * mpi_nProcs
		write(filename, format) eig_out_dir//'ham.',kpt_idx
		!
		open(unit=mpi_unit, file=filename, form='formatted', action='write', access='stream', status='replace')
			do m =1, size(H_k,2)
				do n =1, size(H_k,1)
					write(mpi_unit,'(i2,a,i2,a,e16.7,a,e16.7)') n," ",m," ",real(H_k(n,m),dp)," ",aimag(H_k(n,m))
				end do
			end do
		close(mpi_unit)
		!
		return
	end subroutine


	subroutine write_en_global(kpt_latt)
		real(dp),	intent(in)			::	kpt_latt(:,:)
		real(dp),	allocatable			::	ek_bands(:)
		integer							::	qi_idx, band, x, mpi_unit
		!
		allocate(	ek_bands(num_bands)		)
		!
		mpi_unit	= 200 + mpi_id + 7 * mpi_nProcs
		!
		open(unit=mpi_unit, file=out_dir//'eBands.dat', form='formatted', action='write', access='stream', status='replace')
		write(mpi_unit,*)	'# energies interpolated by MEPinterp program, written '//cTIME(time())
		write(mpi_unit,*)	'# first 3 columns give the relative k-mesh, 4th column are the enegies'
		write(mpi_unit,*)	'# Kpt_idx  K_x (frac)       K_y (frac)        K_z (frac)       Energy (Hartree)  '
		do qi_idx = 1, size(kpt_latt,2)
			call read_en_binary(qi_idx, ek_bands)
			!
			do band = 1, size(ek_bands,1)
				write(mpi_unit,'(i7,a)',advance="no")	qi_idx,'	'
				write(mpi_unit,'(200(f16.8))',advance="no")	(kpt_latt(x,qi_idx), x=1,size(kpt_latt,1)	)
				write(mpi_unit, '(a,f18.8)')			'	',ek_bands(band)
			end do
			!		
		end do
		close(220)
		write(*,'(a,i3,a)')		"[#",mpi_id,"; write_en_global]: success!"
		!
		return
	end subroutine



	subroutine write_velo(kpt_idx, V_ka)
		integer,			intent(in)		::	kpt_idx
		complex(dp),		intent(in)		::	V_ka(:,:,:)
		integer								::	mpi_unit, x, n, m 
		character(len=8)					::	fname
		character(len=120)					::	info_string
		character(len=40)					::	fpath
		!
		!	SPECIFY FILE
		fname	=	'velo_Vka'
		write(info_string,*)	"#interpolated velocities"
		write(fpath,format) 	velo_out_dir //fname//".", kpt_idx 
		write(*,*) fpath
		!
		!	OPEN FILE
		mpi_unit	=	300 + mpi_id + 8 * mpi_nProcs
		open(unit=mpi_unit, file=trim(fpath), form='formatted', action='write', access='stream', status='replace')
		write(mpi_unit,*)	trim(info_string)
		write(mpi_unit,*)	'# x |	 n		|	m 	|	real(V^x_nm)  |		imag(Vx_nm)'
		write(mpi_unit,*)	'#-----------------------------------------------------------------------'
		!
		!	WRITE FILE
		do x = 1, 3
			do m = 1, size(V_ka,3)
				do n = 1, size(V_ka,2)
					write(mpi_unit,'(a,i1,a,i3,a,i3,a,f16.8,a,f16.8)')	"	",x,"		",n,"		",m,"		",	&
																				real(V_ka(x,n,m),dp),"	",aimag(V_ka(x,n,m))
				end do 
			end do
			write(mpi_unit,*)	"#----------------------------------------------------------------------"
		end do
		!
		!	CLOSE FILE
		close(mpi_unit)
		write(*,'(a,i3,a,a,a,i8)')	"[#",mpi_id,";write_velo]: wrote ",fname," at kpt #",kpt_idx
		!
		return
	end subroutine




	subroutine write_hw_list(n_hw, 	hw_min,	hw_max	)
		integer,			intent(in)		::	n_hw
		real(dp),			intent(in)		::	hw_min, hw_max
		integer								::	hw_idx, hw_unit
		real(dp)							::	hw_val, delta_hw
		character(len=20)					::	fname
		character(len=100)					::	info_string
		!
		fname 			=	out_dir//'hw_lst.txt'
		info_string		=	'# List of Laser frequencies (in eV) used for calc. opt. responses, written at '//cTIME(time())
		hw_unit			=	399
		delta_hw		=	(	hw_max -	hw_min	)	/	real(N_hw -1,dp)
		!
		!
		open(unit=hw_unit,	file=trim(fname), form='formatted', action='write', access='stream', status='replace')
		write(hw_unit,*) info_string
		!
		do hw_idx = 1, n_hw
			hw_val		=	hw_min 	+	real(hw_idx-1,dp)	* delta_hw
			!
			write(hw_unit,'(i7,a,f16.8)')	hw_idx,"	", hw_val * aUtoEv
		end do
		close(hw_unit)
		!
		write(*,'(a,i7,a)')		"[write_hw_list]: wrote ",n_hw," entries to ",fname
		!
		return
	end subroutine
















































!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!			WRITE RESPONSE TENSORS:
!--------------------------------------------------------------------------------------------------------------------------------		
!
	subroutine write_mep_bands(n_ki_glob, mep_bands)
		integer,					intent(in)		::	n_ki_glob
		real(dp),	allocatable, 	intent(inout)	::	mep_bands(:,:,:)
		real(dp),	allocatable						::	mep_sum(:,:)
		character(len=20)							::	fname
		character(len=100)							::	info_string
		character(len=3)							::	id_string
		integer										::	n0, row
		logical										::	verbose
		!

		!
		if(allocated(mep_bands)) then
			verbose	= .true.
			allocate(	mep_sum(3,3)	)
			mep_sum		=	0.0_dp
			id_string	=	'mep'
			!
			write(*,*)	"[write_mep_bands]:	"
			do n0 = 1, size(mep_bands,3)
				!	write to file
				write(info_string,'(a,i3,a,i8,a)')		'# mep contribution of band '	,	n0, 	'(n_kpt=',n_ki_glob,')'
				write(fname, format)		 			'mep_band.'						,	n0
				call 	write_tens_file(	mep_out_dir, fname, mep_bands(:,:,n0), info_string, 	id_string,verbose)
				!
				!	print to std output
				write(*,*)	"band n=",n0
				do row = 1, 3
					write(*,*)	mep_bands(row,:,n0)
				end do
				!
				!	add to sum over bands
				mep_sum	=	mep_sum	+	mep_bands(:,:,n0)
			end do
			!
			!	print sum to std output
			write(*,*)	"band sum"
			do row = 1, 3
				write(*,*)	mep_sum(row,:)
			end do
			write(*,*)	"--------------------------"
			write(*,*)	"*"
			write(*,*)	"*"
			!	write sum to file
			write(info_string,*)		'# mep sum over band contributions (this should be the equivalent to mep_tens.dat)'
			fname		='mep_band.sum'
			call 	write_tens_file(	mep_out_dir,	trim(fname),	mep_sum,	info_string,	id_string,verbose)
		end if
		!
		return
	end subroutine


	subroutine write_mep_tensors(n_ki_glob, mep_ic, mep_lc, mep_cs)
		integer,					intent(in)		::	n_ki_glob
		real(dp),	allocatable,	intent(in)		::	mep_ic(:,:), mep_lc(:,:), mep_cs(:,:)
		real(dp),	allocatable						::	mep_tens(:,:)
		character(len=12)							::	fname
		character(len=100)							::	info_string
		character(len=3)							::	id_string
		integer										::	row
		logical										::	verbose
		!
		id_string	=	'mep'
		verbose		=	.true.
		!-------------------------------------total MEP tensor------------------------------
		if(allocated(mep_cs) .and. allocated(mep_lc) .and. allocated(mep_ic) ) then
			fname		= 	'mep_tens.dat'
			info_string	= 	'# total mep tensor (mep_tot= mep_ic+mep_lc+mep_cs)  (in atomic units - dimensionless), written '//cTIME(time())
			!
			allocate(	mep_tens(3,3)	)
			mep_tens	= 	mep_ic +	mep_lc	+	mep_cs
			call	write_tens_file(mep_out_dir,	fname,	mep_tens, info_string,	id_string, verbose )
			write(*,*)	"[write_mep_tensors]: mep_tens:"
			do row = 1, 3
				write(*,*)	mep_tens(row,:)
			end do
		end if
		!
		!-------------------------------------itinerant contribution MEP tensor-------------
		if(allocated(mep_ic)) then
			fname		= 	'mep_ic.dat'
			info_string	= 	'# itinerant contribution of mep tensor (in atomic units - dimensionless), written '//cTIME(time())
			!
			call	write_tens_file(mep_out_dir,	fname,	mep_ic,	info_string,	id_string, verbose )
			!
			write(*,*)	"[write_mep_tensors]: mep_ic:"
			do row = 1, 3
				write(*,*)	mep_ic(row,:)
			end do
		end if
		!
		!
		!-------------------------------------local contribution MEP tensor-----------------
		if(allocated(mep_lc)) then
			fname		= 	'mep_lc.dat'
			info_string	= 	'# local contribution of mep tensor  (in atomic units - dimensionless), written '//cTIME(time())
			!
			call	write_tens_file(mep_out_dir,	fname,	mep_lc,	info_string,	id_string, verbose )
			!
			write(*,*)	"[write_mep_tensors]: mep_lc:"
			do row = 1, 3
				write(*,*)	mep_lc(row,:)
			end do	
		end if
		!
		!
		!-------------------------------------Chern-Simons term MEP tensor------------------
		if(allocated(mep_cs)) then
			fname		= 	'mep_cs.dat'
			info_string	= 	'# Chern-Simons term of mep tensor  (in atomic units - dimensionless), written '//cTIME(time())
			!
			call	write_tens_file(mep_out_dir,	fname,	mep_cs,	info_string,	id_string, verbose )
			!
			write(*,*)	"[write_mep_tensors]: mep_cs:"
			do row = 1, 3
				write(*,*)	mep_cs(row,:)
			end do
		end if
		!-----------------------------------------------------------------------------------
		write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; write_mep_tensors]: calculated MEP tensor on ",n_ki_glob," kpts"
		!
		!
		return
	end subroutine


	subroutine write_kubo_mep_tensors(n_ki_glob, kubo_mep_ic, kubo_mep_lc, kubo_mep_cs)
		integer,					intent(in)		::	n_ki_glob
		real(dp),		allocatable,		intent(in)		::	kubo_mep_ic(:,:), kubo_mep_lc(:,:), kubo_mep_cs(:,:)
		real(dp),		allocatable							::	kubo_mep_tens(:,:)
		character(len=20)									::	fname
		character(len=90)									::	info_string
		character(len=7)									::	id_string
		logical												::	allo_ic, allo_lc, allo_cs, verbose
		!
		id_string	=	'KUBOmep'
		allo_ic		=	.false.
		allo_ic		=	.false.
		allo_ic		=	.false.		
		verbose		=	.true.
		!-------------------------------------itinerant contribution MEP tensor-------------
		if(allocated(kubo_mep_ic)) then
			fname		= 	'kubo_mep_ic.dat'
			info_string	= 	'# itinerant contribution of KUBO mep tensor (with fermi-dirac statistic), written '//cTIME(time())
			call	write_tens_file(mep_out_dir,	fname,	kubo_mep_ic,	info_string,	id_string, verbose )
			allo_ic	=	.true.
		end if
		!
		!
		!-------------------------------------local contribution MEP tensor-----------------
		if(allocated(kubo_mep_lc)) then
			fname		= 	'kubo_mep_lc.dat'
			info_string	= 	'# local contribution of KUBO mep tensor (with fermi-dirac statistic), written '//cTIME(time())
			call	write_tens_file(mep_out_dir,	fname,	kubo_mep_lc,	info_string,	id_string, verbose )
			allo_lc	=	.true.
		end if
		!
		!
		!-------------------------------------Chern-Simons term MEP tensor------------------
		if(allocated(kubo_mep_cs)) then
			fname		= 	'kubo_mep_cs.dat'
			info_string	= 	'# Chern-Simons term of KUBO mep tensor (with fermi-dirac statistic), written '//cTIME(time())
			call	write_tens_file(mep_out_dir,	fname,	kubo_mep_cs,	info_string,	id_string, verbose )
			allo_cs	=	.true.
		end if
		!
		!
		!-------------------------------------total MEP tensor------------------------------
		if(allo_ic .and. allo_lc .and. allo_cs) then
			fname		= 	'kubo_mep_tens.dat'
			info_string	= 	'# total KUBO mep tensor (with fermi-dirac statistic) (mep_tot= mep_ic+mep_lc+mep_cs), written '//cTIME(time())
			allocate(kubo_mep_tens(3,3))
			kubo_mep_tens	= 	kubo_mep_ic +	kubo_mep_lc	+	kubo_mep_cs
			call	write_tens_file(mep_out_dir,	fname,	kubo_mep_tens, info_string,	id_string, verbose )
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; core_worker]: wrote  KUBO MEP tensor on ",n_ki_glob," kpts"
		end if
		!
		!
		return
	end subroutine



	subroutine write_ahc_tensor(n_ki_glob, ahc_tens, velo_ahc_tens, ohc_tens)
		integer,					intent(in)		::	n_ki_glob
		real(dp),	allocatable,	intent(in)			::	ahc_tens(:,:)
		complex(dp),allocatable,	intent(in)			::	velo_ahc_tens(:,:,:), ohc_tens(:,:,:)
		character(len=22)								::	fname
		character(len=100)								::	info_string
		character(len=7)								::	id_string
		integer											::	row, hw_idx, n_hw
		logical											::	verbose
		!
		!
		if(allocated(ahc_tens)) then
			verbose		=	.true.
			fname		=	'ahc_tens.dat'
			info_string	=	'# anomalous Hall conductivity tensor, written '//cTIME(time())
			id_string	=	'ahc'
			call	write_tens_file(ahc_out_dir,	trim(fname),	ahc_tens,	info_string,	id_string, verbose)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; write_ahc_tensor]: wrote AHC tensor (curvature) on ",n_ki_glob," kpts:"
			do row = 1, 3
				write(*,*)	ahc_tens(row,:)
			end do
		end if
		!
		!
		n_hw		=	size(velo_ahc_tens,3)
		!
		verbose	=	.false.
		do hw_idx 	= 1, n_hw
			if(allocated(velo_ahc_tens)) then
				id_string	=	'ahcVELO'
				info_string	=	'# OHC via Wanxiang formula (Im(v**2)/(dE**2 - (hw+eta)**2 )) '//cTIME(time())
				write(fname,format)		'ahc_velo.hw',hw_idx
				
				!
				call	write_tens_file(ahc_out_dir,	fname,	velo_ahc_tens(:,:,hw_idx),	info_string,	id_string, verbose)	
			end if
			!
			if(allocated(ohc_tens)) then
				id_string	=	'ohcVELO'
				info_string	=	'# optical Hall conductivity tensor (via velocities, wann guide: eq.12.5) '//cTIME(time())
				write(fname,format)		'ohc_kubo.hw',hw_idx
				!
				call	write_tens_file(ahc_out_dir,	fname,	ohc_tens(:,:,hw_idx),	info_string,	id_string, verbose)
				!
			
			end if
		end do
		!
		!	DEBUG
		if(allocated(velo_ahc_tens) .and. allocated(ohc_tens)) then
			if(size(velo_ahc_tens,3) /=	size(ohc_tens,3)	) then
				write(*,*)	"[write_ahc_tensor]:	ERROR different hw_lst length detected "
				stop "ahc hw lst inconsistency"
			else
				write(*,'(a,i3,a,i7,a,i8,a)')		"[#",mpi_id,"; write_ahc_tensor]: wrote OHC (Wanxiang) tensor at ",&
																n_hw," frequencies on ",n_ki_glob," kpts"
				write(*,'(a,i3,a,i7,a,i8,a)')		"[#",mpi_id,"; write_ahc_tensor]: wrote OHC (w90: velo) tensor at ",&
																n_hw," frequencies on ",n_ki_glob," kpts"
			end if
		end if
		!
		return
	end subroutine


	subroutine write_opt_tensors(n_ki_glob, s_symm, a_symm, photo_2nd)
		integer,						intent(in)			::	n_ki_glob
		complex(dp),	allocatable,	intent(in)			::	s_symm(:,:,:),	a_symm(:,:,:)
		real(dp),		allocatable,	intent(in)			::	photo_2nd(:,:,:,:)
		character(len=22)									::	fname
		character(len=150)									::	info_string
		character(len=4)									::	id_string
		integer												::	n_hw, hw_idx
		logical												::	verbose
		!
		n_hw 	=	0
		if( 	allocated(photo_2nd		))	n_hw	=	size(	photo_2nd	,	3)
		if( 	allocated(a_symm		)) 	n_hw	=	size(	a_symm		,	3)
		if( 	allocated(s_symm		)) 	n_hw	=	size(	s_symm		,	3)
		!
		if (n_hw > 0) then
			verbose	=	.false.
			do hw_idx	=	 1, n_hw	 
				if(allocated(s_symm)) then
					write(fname,format)		'opt_Ssymm.hw',hw_idx
					info_string	=	'# symmetric optical conductivity tensor, written '//cTIME(time())
					id_string	=	'optS'
					call 	write_tens_file(opt_out_dir,	fname,	s_symm(:,:,hw_idx),	info_string,	id_string, verbose)
				end if
				!
				if(allocated(a_symm)) then
					write(fname,format)		'opt_Asymm.hw',hw_idx
					info_string	=	'# asymmetric optical conductivity tensor, written '//cTIME(time())
					id_string	=	'optA'
					call 	write_tens_file(opt_out_dir,	fname,	a_symm(:,:,hw_idx),	info_string,	id_string, verbose)
				end if
				!
				if(allocated(photo_2nd))	then
					write(fname,format)		'2nd_photo.hw',hw_idx
					write(info_string,*)	'2nd order photcond: J^c_photo = rho^c_ab . E^*_a E_b; Phi_ab= ', 	&
													phi_laser,																&
											' (PRB 97, 241118(R) (2018))'
					id_string	=	"2phC"
					call	write_tens_file(opt_out_dir,	fname,	photo_2nd(:,:,:,hw_idx),	info_string,	id_string, verbose)
				end if
			end do
			!
			!
			write(*,'(a,i3,a,i7,a,i7)')	"[#",mpi_id,";write_opt_tensors]:	wrote optical responses at ",	&
									n_hw," laser frequencies on ",n_ki_glob," kpts"
		end if
		!	
		!
		return
	end subroutine	


	subroutine write_gyro_tensors(n_ki_glob, C_tens, D_tens, Dw_tens)
		integer,					intent(in)		::	n_ki_glob
		complex(dp),	allocatable, intent(in)				::	C_tens(:,:),	D_tens(:,:), Dw_tens(:,:,:)
		character(len=22)									::	fname
		character(len=100)									::	info_string
		character(len=6)									::	id_string
		integer												::	n_hw, hw_idx
		logical												::	verbose
		!
		verbose	=	.true.
		if(allocated(C_tens)) then
			fname		=	'gyro_C.dat'
			info_string	=	'# the C tensor from arXiv:1710.03204v2, written '//cTIME(time())
			id_string	=	'gyroC'
			call	write_tens_file(gyro_out_dir,	fname,	C_tens,	info_string,	id_string, verbose)
		end if
		!
		if(allocated(D_tens)) then
			fname		=	'gyro_D.dat'
			info_string	=	'# the D tensor from arXiv:1710.03204v2, written '//cTIME(time())
			id_string	=	'gyroD'
			call	write_tens_file(gyro_out_dir,	fname,	D_tens,	info_string,	id_string, verbose)
		end if
		!
		if(allocated(Dw_tens)) then
			n_hw		=	size(Dw_tens,3)
			id_string	=	'gyroDw'
			verbose		=	.false.
			!
			do hw_idx = 1, n_hw
				write(fname,format)		'gyro_Dw.hw',hw_idx
				info_string	=	'# the Dw tensor from arXiv:1710.03204v2, written '//cTIME(time())
				call	write_tens_file(gyro_out_dir,	fname,	Dw_tens(:,:,hw_idx),	info_string,	id_string,verbose)
			end do
			!
			write(*,'(a,i3,a,i7,a,i7)')	"[#",mpi_id,";write_gyro_tensors]:	wrote Dw tensor at ",	&
											n_hw,	"laser frequencies on ",n_ki_glob," kpts"		
			!
		end if
		!
		!
		if(allocated(C_tens) .and. allocated(D_tens) .and. allocated(Dw_tens)) then
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; core_worker]: calculated GYRO tensors on ",n_ki_glob," kpts"
		end if
		!
		return
	end subroutine





!private write helpers
	subroutine d2_write_tens_file(dir, fname,tens, info_string, id_string, verbose)
		character(len=*), 	intent(in)		::	dir, fname, info_string, id_string
		real(dp),			intent(in)		::	tens(3,3)
		logical,			intent(in)		::	verbose
		integer								::	row, clm, w_unit		!
		!write result to file
		open(newunit=w_unit, file=dir//fname, form='formatted', 	action='write', access='stream',	status='replace')
			write(w_unit,*)	info_string
			write(w_unit,*)	'begin '//id_string
			do row = 1, 3
				write(w_unit,'(200(f16.8,a))')		(		tens(row,clm), ' ', clm=1,3)
			end do
			write(w_unit,*)	'end '//id_string
		close(w_unit)
		if(verbose)		write(*,'(a,i3,a,a,a,a,a)')	"[#",mpi_id,"; write_",id_string,"_tensor]: file ",&
																fname," written  successfully!"
		!
		!
		return
	end subroutine

	subroutine d3_write_tens_file(dir, fname,tens, info_string, id_string, verbose)
		character(len=*), 	intent(in)		::	dir, fname, info_string, id_string
		real(dp),			intent(in)		::	tens(3,3,3)
		logical,			intent(in)		::	verbose
		integer								::	row, clm,dim, w_unit		!
		character(len=1)					::	dim_str(3)
		!
		dim_str(1)	=	"x"
		dim_str(2)	=	"y"
		dim_str(3)	=	"z"
		!
		!write result to file
		open(newunit=w_unit, file=dir//fname, form='formatted', 	action='write', access='stream',	status='replace')
			write(w_unit,*)	info_string
			write(w_unit,*)	'begin '//id_string
			do dim	=	1,3
				write(w_unit,*)	dim_str(dim)
				do row = 1, 3
					write(w_unit,'(200(f16.8,a))')		(		tens(dim,row,clm), ' ', clm=1,3)
				end do
			end do
			write(w_unit,*)	'end '//id_string
		close(w_unit)
		if(verbose)		write(*,'(a,i3,a,a,a,a,a)')	"[#",mpi_id,"; write_",id_string,"_tensor]: file ",&
																fname," written  successfully!"
		!
		!
		return
	end subroutine



	subroutine z_write_tens_file(dir, fname,tens, info_string, id_string, verbose)
		character(len=*), 	intent(in)		::	dir, fname, info_string, id_string
		complex(dp),		intent(in)		::	tens(3,3)
		logical,			intent(in)		::	verbose
		integer								::	row, clm		!
		!write result to file
		open(unit=255, file=dir//fname, form='formatted', 	action='write', access='stream',	status='replace')
			write(255,*)	info_string
			write(255,*)	'begin '//id_string
			do row = 1, 3
				write(255,'(200(a,f12.6,a,f12.6,a))')		(		' ',real(tens(row,clm),dp), ' ',imag(tens(row,clm)),' ', clm=1,3)
			end do
			write(255,*)	'end '//id_string
		close(255)
		if(verbose)		write(*,'(a,i3,a,a,a,a,a)')	"[#",mpi_id,"; write_",id_string,"_tensor]: file ",&
																fname," written successfully!"
		!
		!
		return
	end subroutine


!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		





















!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!			READ BASIS SET:
!--------------------------------------------------------------------------------------------------------------------------------		
!
!
	logical function read_kptsgen_pl_file(kpt_latt)
		!reads the kpt file generatred by kptsgen.pl
		real(dp),			allocatable,	intent(inout)		::	kpt_latt(:,:)
		character(len=50)										::	filename 
		integer													::	mpi_unit, num_kpts, kpt
		real(dp)												::	raw_kpt(3), tmp
		!
		filename = 'kpts'
		inquire(file=filename, exist= read_kptsgen_pl_file)
		if(read_kptsgen_pl_file ) then
			!
			mpi_unit	= mpi_id + 4*mpi_nProcs
			open(unit=mpi_unit, file=filename, form='formatted', action='read', access='stream', status='old')
			read(mpi_unit,*) num_kpts, raw_kpt(1)
			!
			allocate(	kpt_latt(3,num_kpts)	)
			!
			do kpt = 1, num_kpts
				read(mpi_unit,*) raw_kpt(1:3), tmp
				kpt_latt(1:3,kpt)	= raw_kpt(1:3)
			end do
			close(mpi_unit)
			write(*,'(a,i3,a,i8)')	"[#",mpi_id,"; read_kptsgen_pl_file]: success! num_kpts",num_kpts
		end if 

		!
		return 
	end function


	subroutine read_en_binary(qi_idx, e_bands)
		integer,		intent(in)		::	qi_idx
		real(dp),		intent(out)		::	e_bands(:)
		character(len=24)				::	filepath
		integer							::	mpi_unit
		logical							::	exists
		!
		mpi_unit	= 200 + mpi_id + 5 * mpi_nProcs
		!
		write(filepath, format) raw_dir//'enK.',qi_idx
		inquire(file = filepath, exist=exists)
		if(exists)	then 
			open(unit=mpi_unit, file = filepath, form='unformatted', action='read', access='stream',	status='old'	)
			read(mpi_unit)	e_bands(:)
			close(mpi_unit)
			!write(*,*)	'[read_en_binary]: read ',size(e_bands),' real values from "',trim(filepath),'"'
		else
			e_bands	=	0.0_dp
			write(*,*)	"[read_en_binary]: WARNING could not read ",	filepath,". Does not exist"
		end if
		!
		return
	end subroutine
	

	logical function read_geninterp_kpt_file(seed_name, kpt_latt)
		character(len=*),					intent(in)			::	seed_name
		real(dp),			allocatable,	intent(inout)		::	kpt_latt(:,:)
		character(len=50)										::	filename 
		character(len=500)										::	info_line
		integer													::	mpi_unit, num_kpts, kpt, raw_idx
		real(dp)												::	raw_kpt(3)
		logical													::	found_file, is_frac
		!
		filename	= seed_name//'_geninterp.kpt'
		inquire(file=w90_dir//filename, exist=found_file)
		is_frac	= .false.
		!
		if(	found_file ) then
			mpi_unit	=	mpi_id + 3*mpi_nProcs
			open(unit=mpi_unit, file=w90_dir//filename,	form='formatted', action='read', access='stream', 	status='old'	)
			!
			!read header
			read(mpi_unit,*)
			read(mpi_unit,*)	info_line
			is_frac = 	index(info_line,'frac') /=0
			
			if( is_frac ) then
				read(mpi_unit,*)	num_kpts
				!
				!allocate body container
				allocate(	kpt_latt(3,	num_kpts)	)
				!
				!read body
				do kpt = 1, num_kpts
					read(mpi_unit,*)		raw_idx, raw_kpt(1:3)
					kpt_latt(1:3,raw_idx)	= raw_kpt(1:3)
					if(raw_idx /= kpt)	write(*,'(a,i3,a)') '[#',mpi_id,'; read_geninterp_kpt_file]: WARNING, issues with k-mesh ordering detected'	
				end do
				write(*,'(a,i3,a,i8)')	"[#",mpi_id,"; read_geninterp_kpt_file]: success! num_kpts=",num_kpts
			else
				write(*,'(a,i3,a)')	"[#",mpi_id,"; read_geninterp_kpt_file]: file is not given in fractional coordinates (will not use file)"
			end if
			!
			close(mpi_unit)
		end if
		!
		read_geninterp_kpt_file = found_file .and. is_frac
		return
	end function



	subroutine read_tb_file(seed_name, R_vect, tHopp, rHopp)
		character(len=*)									::	seed_name
		real(dp),		intent(out), allocatable			::	R_vect(:,:)
		complex(dp),	intent(out), allocatable			:: 	tHopp(:,:,:), rHopp(:,:,:,:)
		integer												::	mpi_unit, &
																f_nwfs, f_nSC, readLines, line, &
																cell, n
		integer												::	int3(3), index(2)
		real(dp)											::	unit_cell(3,3), compl1(2),compl3(6), rTest(3), real3(3)
		!
		mpi_unit	= mpi_id
		open(unit=mpi_unit, file=seed_name//'_tb.dat',form='formatted', status='old', action='read')
		!read unit cell
		read(mpi_unit,*)
		read(mpi_unit,*)	real3
		unit_cell(1,1:3)	= real3
		read(mpi_unit,*)	 real3
		unit_cell(2,1:3) 	= real3
		read(mpi_unit,*)	 real3
		unit_cell(3,1:3)	= real3
		unit_cell			= unit_cell / aUtoAngstrm
		!
		!overwrite the lattice read in from input
		a_latt	= unit_cell
		!read basis size
		read(mpi_unit,*)	f_nwfs
		read(mpi_unit,*)	f_nSC
		!allocate targets
		allocate(	R_vect(	3,						f_nSC		)	)
		allocate(	rHopp(	3,	f_nwfs, f_nwfs,		f_nSC		)	)
		allocate(	tHopp(		f_nwfs,	f_nwfs,		f_nSC		)	)
		!dummy read the degeneracy list
		readLines = ceiling( real(f_nSC,dp) / real(15,dp)		)
		do line = 1, readLines
			read(mpi_unit,*)
		end do
		!
		!read tHopp (eV)
		do cell = 1, f_nSC
			read(mpi_unit,*)
			!
			!	READ R-VECTOR AS ...
			if(	use_R_float	) then
				!
				!	... FLOAT
				read(mpi_unit,*)	real3(1:3)
			else
				!
				!	...	INTEGER
				read(mpi_unit,*)	int3(1:3)
				real3(:)	=		int3(:)
				
			end if
			R_vect(1:3,cell)	= matmul(unit_cell,	real3(1:3)	)
			!
			!	READ HOPPING
			do n = 1, f_nWfs**2
				read(mpi_unit,*)	index(1:2), compl1(1:2)
				!
				if (abs(compl1(1))	< 1e-7_dp) 	compl1(1)	=	0.0_dp
				if (abs(compl1(2))	< 1e-7_dp) 	compl1(2)	=	0.0_dp
				!
				tHopp(index(1),index(2), cell)	= cmplx(	compl1(1), compl1(2),	dp		)
			end do
		end do
		!
		!read rHopp (ang)
		do cell = 1, f_nSC
			read(mpi_unit,*)
			!
			!	read R_CELL FROM POS OP
			if(	use_R_float	) then
				read(mpi_unit,*)	real3(1:3)
			else
				read(mpi_unit,*)	int3(1:3)
				real3(:)	=		int3(1:3)
			end if
			rTest(:)	= 	matmul(unit_cell,real3(:)	)
			!
			!	CHECK IF SAME AS ORIGINAL R-CELL FROM HOPPING (see above)
			if(		 norm2( rTest(:) - R_vect(:,cell) )	> 1e-8_dp		) then
				write(*,*)	"[read_tb_basis]: WARNING rHopp has different sc order then tHopp"
			end if
			!
			!
			do n = 1, f_nWfs**2
				read(mpi_unit,*)	index(1:2), compl3(1:6)
				!
				rHopp(1,	index(1), index(2), cell)	= cmplx( compl3(1), compl3(2),	dp	)
				rHopp(2,	index(1), index(2), cell)	= cmplx( compl3(3), compl3(4),	dp	)
				rHopp(3,	index(1), index(2), cell)	= cmplx( compl3(5), compl3(6),	dp	)
			end do
		end do 
		close(mpi_unit)
		!
		!convert to a.u.
		tHopp	= 	tHopp / aUtoEv
		rHopp	= 	rHopp / aUtoAngstrm
		!
		!
		write(*,'(a,i3,a)',advance="no")	"[#",mpi_id,";read_tb_file]: success (input interpretation: nWfs="
		write(*,'(i6,a,i6,a)')				f_nwfs, ";	nrpts=",f_nSC,")!"
		return
	end subroutine 




	subroutine read_hr_file(seed_name, R_vect, H_mat)
		character(len=*),	 			intent(in) 		::	seed_name
		real(dp),		allocatable,	intent(inout)	::	R_vect(:,:)
		complex(dp),	allocatable, 	intent(inout)	::	H_mat(:,:,:)
		integer											::	mpi_unit,&
															f_nwfs, f_nSC, to_read, int15(15), int3(3), m,n, idx, sc, wf
		real(dp)										::	real2(2), real3(3)
		integer,		allocatable						::	R_degneracy(:)
		!
		mpi_unit	= 	100 + mpi_id	+ mpi_nProcs

		open(unit=mpi_unit, file=seed_name//'_hr.dat',	form='formatted', action='read', access='sequential', status='old')
		!write(*,'(a,i3,a,i3)')				 		"[#",mpi_id,";read_hr_file]: 	will open _hr file on unit=",mpi_unit
		!read header
		read(mpi_unit,*)
		read(mpi_unit,*)	f_nwfs
		read(mpi_unit,*)	f_nSC
		!
		allocate(	R_degneracy(					f_nSC	)		)
		allocate(	R_vect(			3,				f_nSC	)		)
		allocate(	H_mat(			f_nWfs, f_nWfs,	f_nSC	)		)
		!
		!read wigner seitz degeneracy
		to_read = f_nSC
		idx		= 1
		do while(to_read>0)
			if(to_read >= 15)	then
				read(mpi_unit,*)		int15(1:15)
				R_degneracy(idx:idx+14)	=	int15(1:15)
				idx = idx+15	
				to_read = to_read-15
			else
				read(mpi_unit,*)		int15(1:to_read)
				R_degneracy(idx:(idx+to_read-1))	= int15(1:to_read)
				idx	= idx + to_read
				to_read	= 0
			end if
		end do
		if(	idx-1 /= f_nSC)	then
			if(mpi_id==mpi_root_id)	 then
				write(*,*)	'idx=',idx,' f_nSC=',f_nSC
				write(*,*)	R_degneracy(1:f_nSC) 
			end if


			stop '[read_hr_file]: ERROR issue reading Wigner Seitz degeneracy (idx/=f_nSC) in  _hr.dat file'
		end if
		!
		!read body
		do sc = 1, f_nSC
			do wf  = 1, f_nwfs**2
				!read next line
				if(use_R_float) then
					read(mpi_unit,*)	real3(1:3),	m, n,	real2(1:2)
				else
					read(mpi_unit,*)		int3(1:3), m,n, 	real2(1:2)
					real3	= real(int3,dp)
				end if
					!
				!get Wigner Seitz vector
				if( wf==1  ) then
					R_vect(1:3,sc)	= real3(1:3)
				else if( .not.	is_equal_vect(fp_acc,	R_vect(1:3,sc),	real3(1:3) )		) then
						write(*,*)				 		"[#",mpi_id,";read_hr_file]: 	WARNING unexpected new R_vect=",real3(:)
						write(*,'(a,i4,a,i4)',advance="no")		"	m=",m," n=",n
						write(*,*)								"	wf=",wf," sc=",sc," expected R_vect=",R_vect(:,sc)
				end if
				!
				!fill Hamiltonian
				if( abs(real2(1))< 1e-7_dp)	real2(1)	=	0.0_dp
				if( abs(real2(2))< 1e-7_dp)	real2(2)	=	0.0_dp
				H_mat(m,n, sc)	= cmplx(	real2(1),	real2(2)	, dp	)
			end do
		end do
		close(mpi_unit)
		!
		!
		!convert to a.u
		H_mat	= H_mat / aUtoEv
		!
		write(*,'(a,i3,a,i4,a,i6,a,i6,a,i5)')	"[#",mpi_id,";read_hr_file]: unit=#",mpi_unit," SUCCESS (input interpretation: nWfs=",		&
												f_nwfs, ";	nrpts=",size(R_vect,2),") on unit: "
		return
	end subroutine



	subroutine read_r_file(seed_name, R_vect, r_mat)
		character(len=*),	 			intent(in) 		::	seed_name
		real(dp),						intent(in)		::	R_vect(:,:)
		complex(dp),	allocatable,	intent(inout)	::	r_mat(:,:,:,:)
		integer											::	mpi_unit, &
															f_nwfs, sc, m, n, int3(3), it
		real(dp)										::	real6(6), real3(3)
		!
		mpi_unit	= 100 + mpi_id + 2*mpi_nProcs
		!
		open(unit=mpi_unit, file=seed_name//'_r.dat',	form='formatted', action='read', access='stream', status='old')
		read(mpi_unit,*)
		read(mpi_unit,*)	f_nwfs
		!
		allocate(	r_mat(3, f_nwfs, f_nwfs,	size(R_vect,2)	)			)
		!
		do sc = 1, size(R_vect,2)
			do it= 1, f_nWfs**2
				if(use_R_float) then
					read(mpi_unit,*)		real3(1:3), 		m, n, 	real6(1:6)
				else
					read(mpi_unit,*)		int3(1:3), 		m, n, 	real6(1:6)
					real3(:)	=	real(int3(:),dp)
				end if
				!
				!
				if( .not.	is_equal_vect(fp_acc,	R_vect(1:3,sc),	real3(1:3))	)	then
					write(*,*)	"[read_r_file]: R_vect=",R_vect(1:3,sc)
					write(*,*)	"[read_r_file]:	input_R=",real(int3(1:3),dp)
					stop '[read_r_file]: different R_vect order in _hr.dat and _r.dat file'
				end if
				!
				r_mat(1,m,n,sc)	=	cmplx(	real6(1)	, real6(2)	, dp	)
				r_mat(2,m,n,sc)	=	cmplx(	real6(3)	, real6(4)	, dp	)
				r_mat(3,m,n,sc)	=	cmplx(	real6(5)	, real6(6)	, dp	)
			end do
		end do
		close(mpi_unit)
		!
		!convert angstom to (au)
		r_mat	= r_mat	/ aUtoAngstrm
		!
		write(*,'(a,i3,a,i6,a)')	"[#",mpi_id,";read_r_file]: success (input interpretation: nWfs=",f_nwfs, ";)"				
		!
		return
	end subroutine



!private helpers
		subroutine read_tb_basis(seed_name, R_vect, H_mat, r_mat)		!ead_tb_basis(seed_name, H_real, r_exist, r_mat)
		!	reads the following files:
		!		'_tb.dat' file			(optional)
		!		if not available
		!			'_hr.dat' file 		(mandatory)
		!			' _r.dat' file		(optional) 
		!
		character(len=*),	 			intent(in) 		::	seed_name
		real(dp),		allocatable,	intent(inout)	::	R_vect(:,:)
		complex(dp),	allocatable, 	intent(inout)	::	H_mat(:,:,:), r_mat(:,:,:,:)
		logical											::	tb_exist, hr_exist, r_exist
		!
		r_exist = .false.
		tb_exist= .false.
		hr_exist= .false.
		!
		!	FULL BASIS
		inquire(file=w90_dir//seed_name//'_tb.dat',	exist=tb_exist)
		if(tb_exist) then
			call read_tb_file(w90_dir//seed_name, R_vect, H_mat, r_mat)
			r_exist	= .true.
		else
			!
			!
			!	HAMILTONIAN
			inquire(file=w90_dir//seed_name//'_hr.dat', exist=hr_exist)
			if( .not. hr_exist)		then
				write(*,*)	'[read_tb_basis]: "'//w90_dir//seed_name//'_hr.dat'//'" does not exist'
				stop '[read_tb_basis]:	ERROR could not read the real space Hamiltonian'
			end if
			call read_hr_file(w90_dir//seed_name, R_vect, H_mat)
			!
			!
			!	POSITION
			inquire(file=w90_dir//seed_name//'_r.dat', exist=r_exist)
			if(  r_exist )	call read_r_file(w90_dir//seed_name, R_vect, r_mat)
			!
			!
		end if
		!
		!	GET BASIS SIZE (NUM_BANDS)
		if(	.not. (size(H_mat,1) == size(H_mat,2)	)	)	stop '[read_tb_basis]:	H_mat is not symmetric'
		num_bands=size(H_mat,1)
		!
		!	CHECK BASIS CONSISTENCY
		if( r_exist	)	then
			if(size(r_mat,1)/=	3			)	stop '[read_tb_basis]:	r_mat does not live in 3D'				
			if(size(r_mat,2) /= size(H_mat,1))	stop '[read_tb_basis]:	r_mat and H_mat have different basis in first dim'
			if(size(r_mat,3) /= size(H_mat,2))	stop '[read_tb_basis]:	r_mat and H_mat have different basis in second dim'
			if(size(r_mat,4) /= size(H_mat,3))	stop '[read_tb_basis]:	r_mat and H_mat live on different real-space (cell) grid'
		end if
		!
		!	EXTENDED DEBUGGING	(HERMITICITY)
		if(debug_mode)	call real_space_debuger(R_vect, H_mat, r_mat)
		!		
		return
	end subroutine


	subroutine real_space_debuger(R_vect, H_mat, r_mat)
		real(dp),						intent(in)			::		R_vect(:,:)	
		complex(dp),	allocatable,	intent(in)			::		H_mat(:,:,:), r_mat(:,:,:,:)
		integer												::		sc, n_cells
		logical												::		hermitian
		real(dp)											::		max_err, R_tot(3)
		!
		R_tot(:) =	0.0_dp
		do sc = 1, size(R_vect,2)
			R_tot(:)	=	R_tot(:)	+	R_vect(:,sc)
		end do 
		write(*,*)	'[read_tb_basis/DEBUG-MODE]:	real space cells (input interpretation):'
		write(*,*)	"	cell	|	Rx 	,	Ry 	, 	Rz	"
		write(*,*)	"----------------------------------------------"
		do sc = 1, size(R_vect,2)
			write(*,'(i3,a,f4.1,a,f4.1,a,f4.1)')	sc,'	|		',	R_vect(1,sc),", ",R_vect(2,sc),', ',R_vect(3,sc)
		end do
		write(*,*)	"----------------------------------------------"
		write(*,'(a,f4.1,a,f4.1,a,f4.1)')	"	sum 	|	",R_tot(1),", ",R_tot(2),', ',R_tot(3)
		if(norm2(R_tot) > 1e-3_dp	) then
			write(*,*)	'[read_tb_basis/DEBUG-MODE]:	WARNING if real space cells dont sum zero hermitictiy issues might arise'
		end if




		!write(*,*)	'[read_tb_basis/DEBUG-MODE]:----start debuging real space basis---------'
		n_cells	=	size(R_vect,2)
		!
		hermitian = size(H_mat,1) == size(H_mat,2)
		do sc = 1, n_cells
			!
			!			HAMILTONIAN
			!
			if( .not.	 is_herm_mat(	H_mat(  :,:,sc),max_err)) then
			 	write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space Hamiltonian is not herm in sc= ',sc,'; largest error: ', max_err
			 	hermitian	=	.false.
			 	!stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (Hamiltonian)"
			end if
			!
			!
			if(allocated(r_mat)) then
				!------------------------------------------------------------------------------------------------------------------------------------
				!			X	- POSITION																											!
				!																																	!
				if( .not.	is_herm_mat(	r_mat(1,:,:,sc),max_err)) then 																			!
					write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space x-position op. not herm, largest error: ', max_err,' sc =',sc			!
					hermitian	=	.false.																											!	
					!stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (x-pos op)"															!	
				end if																																!		
				!																																	!
				!			Y	- POSITION																											!
				!																																	!
				if( .not.	is_herm_mat(	r_mat(2,:,:,sc),max_err)) then 																			!
				write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space y-position op. not herm, largest error: ', max_err	,' sc =',sc							!
					hermitian	=	.false.																											!
					!stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (y-pos op)"															!				
				end if
				!
				!			Z	- POSITION
				!
				if( .not.	is_herm_mat(	r_mat(3,:,:,sc),max_err)) then 
					write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space z-position op. not herm, largest error: ', max_err,' sc =',sc	
					hermitian	=	.false.
					!stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (z-pos op)"
				end if
				!------------------------------------------------------------------------------------------------------------------------------------
			end if

		end do
		!
		if(hermitian)	then
			write(*,*)	"[read_tb_basis/DEBUG-MODE]: SUCCESS real space basis is hermitian"
		else
			write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space basis IS NOT hermitian'
		end if
		!
		!
		!write(*,*)	'[read_tb_basis/DEBUG-MODE]:----finished debuging real space basis---------'
		!
		return
	end subroutine


end module file_io
