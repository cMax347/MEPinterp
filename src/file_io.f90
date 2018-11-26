module file_io
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
	use matrix_math,					only:		is_equal_vect, is_herm_mat	
	use constants,						only:		dp, fp_acc, aUtoAngstrm, aUtoEv
	use mpi_comm,						only:		mpi_id, mpi_root_id, mpi_nProcs
	use input_paras,					only:		w90_dir, 							&
													out_dir,							& 
													velo_out_dir,						&
													mep_out_dir,						&
													ahc_out_dir,						&
													opt_out_dir,						&
													gyro_out_dir,						&
													raw_dir, 							&
													a_latt, 							&
													debug_mode,							&
													plot_bands										

	implicit none

	private
	public									::		read_tb_basis,						& 
													read_kptsgen_pl_file,				&
													!-----------------------------------!
													write_en_binary, read_en_binary,	&
													write_en_global,					&
													write_velo,							&
													!-----------------------------------!
													write_mep_bands,					&
													write_mep_tensors,					&
													write_kubo_mep_tensors,				&
													write_ahc_tensor,					&
													write_opt_tensors,					&
													write_gyro_tensors




	interface write_tens_file
		module procedure d_write_tens_file
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
	subroutine write_en_binary(qi_idx, e_bands)
		integer,		intent(in)		::	qi_idx
		real(dp),		intent(in)		::	e_bands(:)
		character(len=24)				::	filename
		integer							::	mpi_unit
		!
		mpi_unit	=	100 + mpi_id + 6 * mpi_nProcs
		!
		write(filename, format) raw_dir//'enK.',qi_idx
		open(unit=mpi_unit,	file = filename, form='unformatted', action='write', access='stream',	status='replace'		)
		write(mpi_unit)	e_bands(:)
		close(mpi_unit) 
		write(*,*)	"wrote ",size(e_bands), " bands to "//filename 
		!
		return
	end subroutine


	subroutine write_geninterp_kpt_file(seed_name, kpt_latt)
		character(len=*),	intent(in)			::	seed_name
		real(dp),			intent(in)			::	kpt_latt(:,:)	
		integer									::	qi_idx, x		


		open(unit=200, file = seed_name//'_geninterp.kpt', form='formatted', action='write',	access='stream', status='replace')
		write(200,*)	'# fractional kpts created by MEPInterp, written '//cTIME(time())
		write(200,*)	'frac'
		write(200,*)	size(kpt_latt,2)


		do qi_idx = 1, size(kpt_latt,2)
			write(200,'(i3,a)',advance="no")	qi_idx, ' '
			write(200,'(*(f16.8))')		(kpt_latt(x,qi_idx), x= 1, size(kpt_latt,1))
		end do
		close(200)
		write(*,'(a,i3,a,a)')		'[#',mpi_id,'; write_geninterp_kpt_file]: wrote ',seed_name//'_geninterp.kpt file'

		return
	end subroutine



	subroutine write_en_global(n_bands, kpt_latt)
		integer,	intent(in)			::	n_bands
		real(dp),	intent(in)			::	kpt_latt(:,:)
		real(dp),	allocatable			::	ek_bands(:)
		integer							::	qi_idx, band, x, mpi_unit
		!
		allocate(	ek_bands(n_bands)		)
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



	subroutine write_velo(kpt_idx, fname, info_string, V_ka)
		integer,			intent(in)		::	kpt_idx
		character(len=*), 	intent(in)		::	fname
		character(len=*), 	intent(in)		::	info_string
		complex(dp),		intent(in)		::	V_ka(:,:,:)
		integer								::	mpi_unit, x, n, m 
		character(len=40)					::	fpath
		!
		write(fpath,format) velo_out_dir //	trim(fname)//".", kpt_idx 
		!
		mpi_unit	=	300 + mpi_id + 8 * mpi_nProcs
		!
		open(unit=mpi_unit, file=fpath, form='formatted', action='write', access='stream', status='replace')
		write(mpi_unit,*)	trim(info_string)
		write(mpi_unit,*)	'# x |	 n		|	m 	|	real(V^x_nm)  |		imag(Vx_nm)'
		write(mpi_unit,*)	'#-----------------------------------------------------------------------'
		!
		do x = 1, 3
			do m = 1, size(V_ka,3)
				do n = 1, size(V_ka,2)
					write(mpi_unit,'(a,i1,a,i3,a,i3,a,f16.8,a,f16.8)')	"	",x,"		",n,"		",m,"		",	&
																				real(V_ka(x,n,m),dp),"	",aimag(V_ka(x,n,m))
				end do 
			end do
			write(mpi_unit,*)	"#----------------------------------------------------------------------"
		end do
		close(mpi_unit)
		!
		write(*,'(a,i3,a,a,a,i8)')	"[#",mpi_id,";write_velo]: wrote ",fname," at kpt #",kpt_idx
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
		!
		write(*,*)	"WARNING! 	write_mep_bands not implemented yet"
		!
		if(allocated(mep_bands)) then


		end if
		!
		return
	end subroutine




	subroutine write_mep_tensors(n_ki_glob, mep_2014, mep_ic, mep_lc, mep_cs)
		integer,					intent(in)		::	n_ki_glob
		real(dp),	allocatable,	intent(in)		::	mep_2014(:,:), mep_ic(:,:), mep_lc(:,:), mep_cs(:,:)
		real(dp),	allocatable						::	mep_tens(:,:)
		character(len=12)							::	fname
		character(len=100)							::	info_string
		character(len=3)							::	id_string
		integer										::	row
		!
		id_string	=	'mep'
		!-------------------------------------itinerant contribution MEP tensor-------------
		if(allocated(mep_2014)) then
			fname		= 	'mep_14.dat'
			info_string	= 	'# 2014 original paper version of mep tensor, written '//cTIME(time())
			!
			call	write_tens_file(mep_out_dir,	fname,	mep_2014,	info_string,	id_string )
		end if
		!
		!
		!-------------------------------------itinerant contribution MEP tensor-------------
		if(allocated(mep_ic)) then
			fname		= 	'mep_ic.dat'
			info_string	= 	'# itinerant contribution of mep tensor (in atomic units - dimensionless), written '//cTIME(time())
			!
			call	write_tens_file(mep_out_dir,	fname,	mep_ic,	info_string,	id_string )
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
			call	write_tens_file(mep_out_dir,	fname,	mep_lc,	info_string,	id_string )
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
			call	write_tens_file(mep_out_dir,	fname,	mep_cs,	info_string,	id_string )
			!
			write(*,*)	"[write_mep_tensors]: mep_cs:"
			do row = 1, 3
				write(*,*)	mep_cs(row,:)
			end do
		end if
		!
		!
		!-------------------------------------total MEP tensor------------------------------
		if(allocated(mep_cs) .and. allocated(mep_lc) .and. allocated(mep_ic) ) then
			fname		= 	'mep_tens.dat'
			info_string	= 	'# total mep tensor (mep_tot= mep_ic+mep_lc+mep_cs)  (in atomic units - dimensionless), written '//cTIME(time())
			!
			allocate(	mep_tens(3,3)	)
			mep_tens	= 	mep_ic +	mep_lc	+	mep_cs
			call	write_tens_file(mep_out_dir,	fname,	mep_tens, info_string,	id_string )
			write(*,*)	"[write_mep_tensors]: mep_tens:"
			do row = 1, 3
				write(*,*)	mep_tens(row,:)
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
		logical												::	allo_ic, allo_lc, allo_cs
		!
		id_string	=	'KUBOmep'
		allo_ic		=	.false.
		allo_ic		=	.false.
		allo_ic		=	.false.		
		!-------------------------------------itinerant contribution MEP tensor-------------
		if(allocated(kubo_mep_ic)) then
			fname		= 	'kubo_mep_ic.dat'
			info_string	= 	'# itinerant contribution of KUBO mep tensor (with fermi-dirac statistic), written '//cTIME(time())
			call	write_tens_file(mep_out_dir,	fname,	kubo_mep_ic,	info_string,	id_string )
			allo_ic	=	.true.
		end if
		!
		!
		!-------------------------------------local contribution MEP tensor-----------------
		if(allocated(kubo_mep_lc)) then
			fname		= 	'kubo_mep_lc.dat'
			info_string	= 	'# local contribution of KUBO mep tensor (with fermi-dirac statistic), written '//cTIME(time())
			call	write_tens_file(mep_out_dir,	fname,	kubo_mep_lc,	info_string,	id_string )
			allo_lc	=	.true.
		end if
		!
		!
		!-------------------------------------Chern-Simons term MEP tensor------------------
		if(allocated(kubo_mep_cs)) then
			fname		= 	'kubo_mep_cs.dat'
			info_string	= 	'# Chern-Simons term of KUBO mep tensor (with fermi-dirac statistic), written '//cTIME(time())
			call	write_tens_file(mep_out_dir,	fname,	kubo_mep_cs,	info_string,	id_string )
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
			call	write_tens_file(mep_out_dir,	fname,	kubo_mep_tens, info_string,	id_string )
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; core_worker]: wrote  KUBO MEP tensor on ",n_ki_glob," kpts"
		end if
		!
		!
		return
	end subroutine


	subroutine write_ahc_tensor(n_ki_glob, ahc_tens, velo_ahc_tens, ohc_tens)
		integer,					intent(in)		::	n_ki_glob
		real(dp),	allocatable,	intent(in)			::	ahc_tens(:,:), velo_ahc_tens(:,:)
		complex(dp),allocatable,	intent(in)			::	ohc_tens(:,:)
		character(len=12)								::	fname
		character(len=70)								::	info_string
		character(len=7)								::	id_string
		!
		!
		if(allocated(ahc_tens)) then
			fname		=	'ahc_tens.dat'
			info_string	=	'# anomalous Hall conductivity tensor, written '//cTIME(time())
			id_string	=	'ahc'
			call	write_tens_file(ahc_out_dir,	fname,	ahc_tens,	info_string,	id_string)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; write_ahc_tensor]: wrote AHC tensor on ",n_ki_glob," kpts"
		end if
		!
		!
		if(allocated(velo_ahc_tens)) then
			fname		=	'ahc_velo.dat'
			info_string	=	'# anomalous Hall conductivity tensor (via velocity kubo), written '//cTIME(time())
			id_string	=	'ahcVELO'
			call	write_tens_file(ahc_out_dir,	fname,	velo_ahc_tens,	info_string,	id_string)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; write_ahc_tensor]: wrote AHC (via velo) tensor on ",n_ki_glob," kpts"
		end if
		!
		!
		if(allocated(ohc_tens)) then
			fname		=	'ohc_kubo.dat'
			info_string	=	'# optical Hall conductivity tensor (wann guide: eq.12.5) '//cTIME(time())
			id_string	=	'ohcVELO'
			call	write_tens_file(ahc_out_dir,	fname,	ohc_tens,	info_string,	id_string)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; write_ahc_tensor]: wrote OHC (via velo) tensor on ",n_ki_glob," kpts"
		end if
		!
		!
		return
	end subroutine


	subroutine write_opt_tensors(n_ki_glob, s_symm, a_symm)
		integer,					intent(in)		::	n_ki_glob
		complex(dp),	allocatable,	intent(in)			::	s_symm(:,:),	a_symm(:,:)
		character(len=13)									::	fname
		character(len=70)									::	info_string
		character(len=4)									::	id_string
		!
		!
		if(allocated(s_symm)) then
			fname		=	'opt_Ssymm.dat'
			info_string	=	'# symmetric optical conductivity tensor, written '//cTIME(time())
			id_string	=	'optS'
			call 	write_tens_file(opt_out_dir,	fname,	s_symm,	info_string,	id_string)
		end if
		!
		if(allocated(a_symm)) then
			fname		=	'opt_Asymm.dat'
			info_string	=	'# asymmetric optical conductivity tensor, written '//cTIME(time())
			id_string	=	'optA'
			call 	write_tens_file(opt_out_dir,	fname,	a_symm,	info_string,	id_string)
			!
			!
			if(allocated(s_symm)) then
				write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; core_worker]: calculated OPT tensor on ",n_ki_glob," kpts"
			end if
		end if
		!
		!
		return
	end subroutine	


	subroutine write_gyro_tensors(n_ki_glob, C_tens, D_tens, Dw_tens)
		integer,					intent(in)		::	n_ki_glob
		complex(dp),	allocatable, intent(in)				::	C_tens(:,:),	D_tens(:,:), Dw_tens(:,:)
		character(len=13)									::	fname
		character(len=70)									::	info_string
		character(len=6)									::	id_string
		!
		if(allocated(C_tens)) then
			fname		=	'gyro_C.dat'
			info_string	=	'# the C tensor from arXiv:1710.03204v2, written '//cTIME(time())
			id_string	=	'gyroC'
			call	write_tens_file(gyro_out_dir,	fname,	C_tens,	info_string,	id_string)
		end if
		!
		if(allocated(D_tens)) then
			fname		=	'gyro_D.dat'
			info_string	=	'# the D tensor from arXiv:1710.03204v2, written '//cTIME(time())
			id_string	=	'gyroD'
			call	write_tens_file(gyro_out_dir,	fname,	D_tens,	info_string,	id_string)
		end if
		!
		if(allocated(Dw_tens)) then
			fname		=	'gyro_Dw.dat'
			info_string	=	'# the Dw tensor from arXiv:1710.03204v2, written '//cTIME(time())
			id_string	=	'gyroDw'
			call	write_tens_file(gyro_out_dir,	fname,	Dw_tens,	info_string,	id_string)
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
	subroutine d_write_tens_file(dir, fname,tens, info_string, id_string)
		character(len=*), 	intent(in)		::	dir, fname, info_string, id_string
		real(dp),			intent(in)		::	tens(3,3)
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
		write(*,'(a,i3,a,a,a,a,a)')	"[#",mpi_id,"; write_",id_string,"_tensor]: file ",fname," written  successfully!"
		!
		!
		return
	end subroutine


	subroutine z_write_tens_file(dir, fname,tens, info_string, id_string)
		character(len=*), 	intent(in)		::	dir, fname, info_string, id_string
		complex(dp),		intent(in)		::	tens(3,3)
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
		write(*,'(a,i3,a,a,a,a,a)')	"[#",mpi_id,"; write_",id_string,"_tensor]: file ",fname," written successfully!"
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
!			READ K-SPACE :
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
			write(*,*)	'[read_en_binary]: read ',size(e_bands),' real values from "',trim(filepath),'"'
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








!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!			READ REAL SPACE :
!--------------------------------------------------------------------------------------------------------------------------------		
!
!






	subroutine read_tb_file(seed_name, R_vect, tHopp, rHopp)
		character(len=*)									::	seed_name
		real(dp),		intent(out), allocatable			::	R_vect(:,:)
		complex(dp),	intent(out), allocatable			:: 	tHopp(:,:,:), rHopp(:,:,:,:)
		integer												::	mpi_unit, &
																f_nwfs, f_nSC, readLines, line, &
																cell, n
		integer												::	cell_rel(3), index(2)
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
			read(mpi_unit,*) cell_rel(1:3)

			R_vect(1:3,cell)	= matmul(unit_cell,	real(cell_rel(1:3),dp)	)
			do n = 1, f_nWfs**2
				read(mpi_unit,*)	index(1:2), compl1(1:2)
				!
				tHopp(index(1),index(2), cell)	= cmplx(	compl1(1), compl1(2),	dp		)
			end do
		end do
		!
		!read rHopp (ang)
		do cell = 1, f_nSC
			read(mpi_unit,*)
			read(mpi_unit,*)	cell_rel(1:3)

			rTest(1:3)	= matmul(unit_cell,	real(cell_rel(1:3),dp)	)
			if(		 norm2( rTest(1:3) - R_vect(1:3,cell) )	> 1e-8_dp		) then
				write(*,*)	"[read_tb_basis]: WARNING rHopp has different sc order then tHopp"
			end if

			do n = 1, f_nWfs**2
				read(mpi_unit,*)	index(1:2), compl3(1:6)

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
		real(dp)										::	real2(2)
		integer,		allocatable						::	R_degneracy(:)
		!
		mpi_unit	= 	100 + mpi_id	+ mpi_nProcs

		open(unit=mpi_unit, file=seed_name//'_hr.dat',	form='formatted', action='read', access='sequential', status='old')
		write(*,'(a,i3,a,i3)')				 		"[#",mpi_id,";read_hr_file]: 	will open _hr file on unit=",mpi_unit
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


			stop 'issue reading Wigner Seitz degeneracy (idx/=f_nSC) in  _hr.dat file'
		end if
		!
		!read body
		idx = 0
		do sc = 1, f_nSC
			do wf  = 1, f_nwfs**2
				!read next line
				read(mpi_unit,*)		int3(1:3), m,n, 	real2(1:2)
				!
				!get Wigner Seitz vector
				if( wf==1 .and. sc==1 ) then
					idx = 1
					R_vect(1:3,idx)	= real(int3(1:3),dp)
				else if( .not.	is_equal_vect(fp_acc,	R_vect(1:3,idx),	real(int3(1:3),dp) )		) then
					idx = idx +1 
					R_vect(1:3,idx)	= real(int3(1:3),dp)
					if(wf /= 1)	then 
						write(*,'(a,i3,a)')				 		"[#",mpi_id,";read_hr_file]: 	WARNING unexpected new R_vect"
						write(*,'(a,i4,a,i4)',advance="no")		"	m=",m," n=",n
						write(*,*)								"	wf=",wf," sc=",sc," raw R_vect=",R_vect(:,idx)
					end if
				end if
				!
				!fill Hamiltonian
				H_mat(m,n, idx)	= cmplx(	real2(1),	real2(2)	, dp	)
			end do
		end do
		close(mpi_unit)
		if(idx /= f_nSC)	stop 'issue reading the body in _hr.dat file'
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
		real(dp)										::	real6(6)
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
				read(mpi_unit,*)		int3(1:3), 		m, n, 	real6(1:6)
				if( .not.	is_equal_vect(fp_acc,	R_vect(1:3,sc),	real(int3(1:3),dp))	)	then
					write(*,*)	"[read_r_file]: R_vect=",R_vect(1:3,sc)
					write(*,*)	"[read_r_file]:	input_R=",real(int3(1:3),dp)
					stop 'different R_vect order in _hr.dat and _r.dat file'
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
			 	stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (Hamiltonian)"
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
					stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (x-pos op)"															!	
				end if																																!		
				!																																	!
				!			Y	- POSITION																											!
				!																																	!
				if( .not.	is_herm_mat(	r_mat(2,:,:,sc),max_err)) then 																			!
				write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space y-position op. not herm, largest error: ', max_err	,' sc =',sc							!
					hermitian	=	.false.																											!
					stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (y-pos op)"															!				
				end if
				!
				!			Z	- POSITION
				!
				if( .not.	is_herm_mat(	r_mat(3,:,:,sc),max_err)) then 
					write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space z-position op. not herm, largest error: ', max_err,' sc =',sc	
					hermitian	=	.false.
					stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (z-pos op)"
				end if
				!------------------------------------------------------------------------------------------------------------------------------------
			end if

		end do
		!
		if(hermitian)	then
			write(*,*)	"[read_tb_basis/DEBUG-MODE]: SUCCESS real space basis is hermitian"
		else
			stop '[read_tb_basis/DEBUG-MODE]: ERROR real space basis IS NOT hermitian'
		end if
		!
		!
		!write(*,*)	'[read_tb_basis/DEBUG-MODE]:----finished debuging real space basis---------'
		!
		return
	end subroutine


end module file_io
