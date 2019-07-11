module file_io
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
	use m_npy
	use endian_swap
	use matrix_math,					only:		is_equal_vect, is_herm_mat	
	use constants,						only:		dp, fp_acc, aUtoAngstrm, aUtoEv
	use mpi_community,					only:		mpi_id, mpi_root_id, mpi_nProcs
	use input_paras,					only:		w90_dir, 							&
													out_dir,							& 
													velo_out_dir,						&
													raw_dir, 							&
													a_latt, 							&
													debug_mode,							&
													use_R_float,						&
													plot_bands									

	implicit none

	private
	public									::		write_ham_binary,					&
													write_eig_binary,					&
													write_en_binary, read_en_binary,	&
													write_en_global,					&
													write_velo

	character(len=64)						::		format='(a,i7.7)'
	


	contains





!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!			WRITE ROUTINES
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
	!~
	!~
	subroutine write_eig_binary(kpt_idx, U_k)
		integer,		intent(in)		::	kpt_idx
		complex(dp),	intent(in)		::	U_k(:,:)
		character(len=24)				::	filename
		integer							::	mpi_unit, n, m
		!
		mpi_unit	=	100 + mpi_id + 15 * mpi_nProcs
		write(filename, format) raw_dir//'eig.',kpt_idx
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
	!~
	!~
	subroutine write_ham_binary(kpt_idx, H_k)
		integer,		intent(in)		::	kpt_idx
		complex(dp),	intent(in)		::	H_k(:,:)
		character(len=24)				::	filename
		integer							::	mpi_unit, n, m
		!
		mpi_unit	=	100 + mpi_id + 15 * mpi_nProcs
		write(filename, format) raw_dir//'ham.',kpt_idx
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
	!~
	!~
	subroutine write_en_global(kpt_latt, num_bands)
		real(dp),	intent(in)			::	kpt_latt(:,:)
		integer,	intent(in)			::	num_bands
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
			do band = 1, size(ek_bands,1)
				write(mpi_unit,'(i7,a)',advance="no")	qi_idx,'	'
				write(mpi_unit,'(200(f16.8))',advance="no")	(kpt_latt(x,qi_idx), x=1,size(kpt_latt,1)	)
				write(mpi_unit, '(a,f18.8)')			'	',ek_bands(band)
			end do
			!		
		end do
		close(220)
		write(*,'(a,i7.7,a)')		"[#",mpi_id,"; write_en_global]: success!"
		!
		return
	end subroutine
	!~
	!~
	subroutine write_velo(kpt_idx, V_ka)
		integer,				 	intent(in)		::	kpt_idx
		complex(dp), allocatable,	intent(in)		::	V_ka(:,:,:)
		integer										::	mpi_unit, x, n, m 
		character(len=8)							::	fname
		character(len=120)							::	info_string
		character(len=40)							::	fpath
		!
		if(allocated(V_ka)) then
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
			write(*,'(a,i7.7,a,a,a,i8)')	"[#",mpi_id,";write_velo]: wrote ",fname," at kpt #",kpt_idx
		else
			write(*,'(a,i7.7,a)')	"[#",mpi_id,";write_velo]: ERROR velocity array not allocated"
		end if
		!
		return
	end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!--------------------------------------------------------------------------------------------------------------------------------		
!			READ ROUTINES
!--------------------------------------------------------------------------------------------------------------------------------		
!
!



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
		write(*,*)	"[read_en_binary]: read ",size(e_bands)," bands for #qi",qi_idx
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
					if(raw_idx /= kpt)	write(*,'(a,i7.7,a)') '[#',mpi_id,&
						'; read_geninterp_kpt_file]: WARNING, issues with k-mesh ordering detected'	
				end do
				write(*,'(a,i7.7,a,i8)')	"[#",mpi_id,"; read_geninterp_kpt_file]: success! num_kpts=",num_kpts
			else
				write(*,'(a,i7.7,a)')	"[#",mpi_id,"; read_geninterp_kpt_file]: file is not given in fractional coordinates (will not use file)"
			end if
			!
			close(mpi_unit)
		end if
		!
		read_geninterp_kpt_file = found_file .and. is_frac
		return
	end function


!~
!~
!~
!~
end module file_io







!!--------------------------------------------------------------------------------------------------------------------------------		
!!--------------------------------------------------------------------------------------------------------------------------------		
!!--------------------------------------------------------------------------------------------------------------------------------		
!!			WRITE RESPONSE TENSORS:	(REPLACED BY DIRECT SAVE_NPY CALL IN CORE WORKER!!)
!!--------------------------------------------------------------------------------------------------------------------------------		
!!
!	subroutine write_mep_bands(n_ki_glob, mep_bands)
!		integer,					intent(in)		::	n_ki_glob
!		real(dp),	allocatable, 	intent(inout)	::	mep_bands(:,:,:)
!		!
!		if(	allocated(mep_bands))	then
!			call	save_npy(out_dir//"mep_bands.npy",	mep_bands	)
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,"; write_mep_tensors]: wrote mep_bands to ",out_dir//"mep_bands.npy"
!		end if
!		!
!		return
!	end subroutine
!	!
!	subroutine write_mep_tensors(n_ki_glob, mep_ic, mep_lc, mep_cs)
!		integer,					intent(in)		::	n_ki_glob
!		real(dp),	allocatable,	intent(in)		::	mep_ic(:,:), mep_lc(:,:), mep_cs(:,:)
!		real(dp)									::	mep_tens(3,3)
!		!
!		if(allocated(mep_ic)) 	call save_npy(	out_dir//'mep_ic.npy', 		mep_ic		)
!		if(allocated(mep_lc)) 	call save_npy(	out_dir//'mep_lc.npy', 		mep_lc		)
!		if(allocated(mep_cs)) 	call save_npy(	out_dir//'mep_cs.npy', 		mep_cs		)
!		!
!		if(allocated(mep_ic)) 	&
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,"; write_mep_tensors]: wrote mep to ",out_dir//'mep_ic.npy'
!		if(allocated(mep_lc))	&
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,"; write_mep_tensors]: wrote mep to ",out_dir//'mep_lc.npy'
!		if(allocated(mep_cs))	&
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,"; write_mep_tensors]: wrote mep to ",out_dir//'mep_cs.npy'
!		!-------------------------------------total MEP tensor------------------------------
!		if(allocated(mep_cs) .and. allocated(mep_lc) .and. allocated(mep_ic) ) then
!			mep_tens	= 	mep_ic +	mep_lc	+	mep_cs
!			call save_npy(	out_dir//'mep_tens.npy',		mep_tens	)
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,"; write_mep_tensors]: wrote mep to ",out_dir//'mep_tens.npy'
!		end if
!		!-----------------------------------------------------------------------------------
!		
!		!
!		!
!		return
!	end subroutine
!	!
!	subroutine write_kubo_mep_tensors(n_ki_glob, kubo_mep_ic, kubo_mep_lc, kubo_mep_cs)
!		integer,							intent(in)		::	n_ki_glob
!		real(dp),		allocatable,		intent(in)		::	kubo_mep_ic(:,:), kubo_mep_lc(:,:), kubo_mep_cs(:,:)
!		real(dp)											::	kubo_mep_tens(3,3)
!		!
!		!-------------------------------------individual contributions-------------
!		if(allocated(kubo_mep_ic)) 	call save_npy(	out_dir//'kubo_mep_ic.npy',	kubo_mep_ic	)
!		if(allocated(kubo_mep_lc)) 	call save_npy(	out_dir//'kubo_mep_lc.npy',	kubo_mep_lc	)
!		if(allocated(kubo_mep_cs)) 	call save_npy(	out_dir//'kubo_mep_cs.npy',	kubo_mep_cs	)
!		!
!		if(allocated(kubo_mep_ic))	&
!			write(*,'(a,i7.7,a,a)')	"[#",mpi_id,"; write_mep_tensors]: wrote kubo_mep to ",out_dir//'kubo_mep_ic.npy'
!		if(allocated(kubo_mep_lc))	&
!			write(*,'(a,i7.7,a,a)')	"[#",mpi_id,"; write_mep_tensors]: wrote kubo_mep to ",out_dir//'kubo_mep_lc.npy'
!		if(allocated(kubo_mep_cs))	&
!			write(*,'(a,i7.7,a,a)')	"[#",mpi_id,"; write_mep_tensors]: wrote kubo_mep to ",out_dir//'kubo_mep_cs.npy'
!		!
!		!-------------------------------------total MEP tensor------------------------------
!		if(			allocated(kubo_mep_ic) 		&	
!			.and. 	allocated(kubo_mep_lc)		&
!			.and. 	allocated(kubo_mep_cs)		&
!		) then
!			kubo_mep_tens	= 	kubo_mep_ic +	kubo_mep_lc	+	kubo_mep_cs
!			call save_npy(	out_dir//'kubo_mep_tens.npy',	kubo_mep_tens	)
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,"; write_mep_tensors]: wrote kubo_mep to ",out_dir//'kubo_mep_tens.npy'
!		end if
!		!
!		return
!	end subroutine
!	!
!	subroutine write_ahc_tensor(n_ki_glob, ahc_tens, velo_ahc_tens, ohc_tens)
!		integer,					intent(in)		::	n_ki_glob
!		real(dp),	allocatable,	intent(in)			::	ahc_tens(:,:)
!		complex(dp),allocatable,	intent(in)			::	velo_ahc_tens(:,:,:), ohc_tens(:,:,:)
!		!
!		!
!		if(allocated(ahc_tens)) 	call 	save_npy(	out_dir//'ahc_tens.npy',	ahc_tens		)
!		if(allocated(velo_ahc_tens))call	save_npy(	out_dir//'ahcVELO.npy',		velo_ahc_tens	)
!		if(allocated(ohc_tens))		call	save_npy(	out_dir//'ohcVELO.npy',		ohc_tens		)
!		!
!		!
!		if(allocated(ahc_tens)) 	&
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,"; write_ahc_tensor]: wrote AHC/OHC tensor to ",out_dir//'ahc_tens.npy'
!		if(allocated(velo_ahc_tens))&
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,"; write_ahc_tensor]: wrote AHC/OHC tensor to ",out_dir//'ahcVELO.npy'
!		if(allocated(ohc_tens))		&
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,"; write_ahc_tensor]: wrote AHC/OHC tensor to ",out_dir//'ohcVELO.npy'
!		!
!		return
!	end subroutine
!	!
!	subroutine write_opt_tensors(n_ki_glob, s_symm, a_symm)
!		integer,						intent(in)			::	n_ki_glob
!		complex(dp),	allocatable,	intent(in)			::	s_symm(:,:,:),	a_symm(:,:,:)
!		
!		!
!		if(allocated(s_symm)) 		call save_npy(	out_dir//"opt_Ssymm.npy", 	s_symm		)
!		if(allocated(a_symm))		call save_npy(	out_dir//"opt_Asymm.npy", 	a_symm		)
!		
!		!
!		if(allocated(s_symm))	&
!			write(*,'(a,i7.7,a,a)')	"[#",mpi_id,";write_opt_tensors]:	wrote symm opt tens to ",	out_dir//"opt_Ssymm.npy"
!		if(allocated(a_symm))	&
!			write(*,'(a,i7.7,a,a)')	"[#",mpi_id,";write_opt_tensors]:	wrote symm opt tens to ",	out_dir//"opt_Asymm.npy"
!		return
!	end subroutine	
!	!
!	subroutine write_photoC_tensors(n_ki_glob, photo_2nd)
!		integer,						intent(in)			::	n_ki_glob
!		complex(dp),	allocatable,	intent(in)			::	photo_2nd(:,:,:,:,:,:)
!		!
!		if(allocated(photo_2nd))	call save_npy(	out_dir//'photoC_2nd.npy', 	photo_2nd	)
!		if(allocated(photo_2nd))	&
!			write(*,'(a,i7.7,a,a)')	"[#",mpi_id,";write_opt_tensors]:	wrote 2nd opt tens to ",	out_dir//'photoC_2nd.npy'
!		!
!		return
!	end subroutine
!	!
!	subroutine write_gyro_tensors(n_ki_glob, C_tens, D_tens, Dw_tens)
!		integer,					intent(in)		::	n_ki_glob
!		complex(dp),	allocatable, intent(in)		::	C_tens(:,:,:),	D_tens(:,:,:), Dw_tens(:,:,:,:)
!		!
!		if(allocated(C_tens)) 		call save_npy(	out_dir//"gyro_C.npy", 	C_tens		)
!		if(allocated(D_tens)) 		call save_npy(	out_dir//"gyro_D.npy", 	D_tens		)
!		if(allocated(Dw_tens))		call save_npy(	out_dir//"gyro_Dw.npy", Dw_tens		)
!		!
!		if(allocated(C_tens))	&
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,";write_gyro_tensors]: wrote GYRO-C tensors to ",out_dir//"gyro_C.npy"
!		if(allocated(D_tens))	&
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,";write_gyro_tensors]: wrote GYRO-D tensors to ",out_dir//"gyro_D.npy"
!		if(allocated(Dw_tens))	&	
!			write(*,'(a,i7.7,a,a)')		"[#",mpi_id,";write_gyro_tensors]: wrote GYR-DW tensors to ",out_dir//"gyro_Dw.npy"
!		!
!		return
!	end subroutine
!!--------------------------------------------------------------------------------------------------------------------------------		
!!--------------------------------------------------------------------------------------------------------------------------------		
!!--------------------------------------------------------------------------------------------------------------------------------		
