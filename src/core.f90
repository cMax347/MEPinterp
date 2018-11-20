module core
	!this module uses a semiclassic approach to calculate the first ordrer correction
	!	to the polariztion induced by a perturbive magnetic field
	! 	see Niu PRL 112, 166601 (2014)
	!use omp_lib
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
#ifdef USE_MPI
	use mpi
#endif
	use matrix_math,	only:	my_Levi_Civita
	use constants,		only:	dp
	use mpi_comm,		only:	mpi_root_id, mpi_id, mpi_nProcs, ierr,			&
								mpi_ki_selector,								&
								mpi_reduce_tens
	use statistics,		only:	fd_count_el
	use input_paras,	only:	use_mpi,										&
								do_gauge_trafo,									&
								do_mep,											&
								do_write_mep_bands,								&
								do_kubo,										&
								do_ahc,											&
								do_opt,											&
								do_gyro,										&
								debug_mode,										&
								a_latt, 										&
								valence_bands, 									&
								seed_name,										&
								kubo_tol,										&
								hw, eFermi, T_kelvin, i_eta_smr, 				&
								unit_vol
	!
	use k_space,		only:	get_recip_latt, get_mp_grid, 					&
								kspace_allocator,								&
								get_rel_kpt,									&
								normalize_k_int
	!								
	use wann_interp,	only:	get_wann_interp
	!
	use mep_niu,		only:	mep_niu_CS,	mep_niu_IC, mep_niu_LC
	use mep_niu2014,	only:	mep_niu2014_full
	use kubo,			only:	kubo_ahc_tens, velo_ahc_tens, kubo_opt_tens
	use kubo_mep,		only:	kubo_mep_CS, kubo_mep_LC, kubo_mep_IC
	use gyro,			only:	get_gyro_C, get_gyro_D, get_gyro_Dw
	!
	use file_io,		only:	read_tb_basis,									&
								write_mep_bands,								&
								write_mep_tensors,								&
								write_kubo_mep_tensors,							&
								write_ahc_tensor,								&
								write_opt_tensors,								&
								write_gyro_tensors
	


	implicit none



	private
	public ::			core_worker
	!
	save
	integer									::		num_kpts
contains



!public
	subroutine	core_worker()
		!	interpolate the linear magnetoelectric response tensor alpha
		!
		!	lc	: local current contribution
		!	ic	: itinerant current contribution
		!	cs	: chern - simons term	
		!
		real(dp)							::	N_el_k,  max_n_el, min_n_el,			&
												sum_N_el_loc, 							&
												kpt(3),	recip_latt(3,3)
												!local sum targets:
		real(dp),		allocatable			::	mep_2014_loc(		:,:),				&
												kubo_ahc_loc(		:,:),				&
												velo_ahc_loc(		:,:),				&
												kubo_mep_ic_loc(	:,:),				&
												kubo_mep_lc_loc(	:,:),				&
												kubo_mep_cs_loc(	:,:)
												!
		complex(dp),	allocatable			::	tempS(:,:), 	tempA(:,:),				&
												kubo_opt_s_loc(		:,:),				&
												kubo_opt_a_loc(		:,:),				&
												gyro_D_loc(			:,:),				&
												gyro_Dw_loc(		:,:),				&																								
												gyro_C_loc(			:,:)																								
												!
		integer								::	kix, kiy, kiz, ki, 						&
												mp_grid(3), n_ki_loc, 					&
												ic_skipped, lc_skipped
		complex(dp),	allocatable			::	H_tb(:,:,:), r_tb(:,:,:,:), 			&
												A_ka(:,:,:), Om_kab(:,:,:,:),			&
												V_ka(:,:,:)
		real(dp),		allocatable			::	en_k(:), R_vect(:,:),					& 
												mep_tens_ic_loc(	:,:,:),	 			&
												mep_tens_lc_loc(	:,:,:),				&
												mep_tens_cs_loc(	:,:,:)
		!----------------------------------------------------------------------------------------------------------------------------------
		!	allocate
		!----------------------------------------------------------------------------------------------------------------------------------
		call allo_core_loc_arrays(	N_el_k, sum_N_el_loc, max_n_el, min_n_el,										&
									mep_2014_loc,																	&
									mep_tens_ic_loc, mep_tens_lc_loc, mep_tens_cs_loc,								&	
									kubo_mep_ic_loc, kubo_mep_lc_loc, kubo_mep_cs_loc,								&									
									kubo_ahc_loc, velo_ahc_loc,														&
									tempS, tempA, kubo_opt_s_loc, kubo_opt_a_loc,									&
									gyro_C_loc, gyro_D_loc, gyro_Dw_loc												&
						)
		!
		!----------------------------------------------------------------------------------------------------------------------------------
		!	get 	k-space
		!----------------------------------------------------------------------------------------------------------------------------------
		mp_grid				=	get_mp_grid()
		n_ki_loc			= 	0
		num_kpts			= 	mp_grid(1)*mp_grid(2)*mp_grid(3)
		recip_latt			=	get_recip_latt()
		!
		!----------------------------------------------------------------------------------------------------------------------------------
		!	get		TB	BASIS
		!----------------------------------------------------------------------------------------------------------------------------------
		if(use_mpi) call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		if(mpi_id==mpi_root_id)	then
			write(*,*)	"*"
			write(*,*)	"----------------------GET REAL SPACE BASIS---------------------------"
		end if
		call read_tb_basis(			seed_name, R_vect,		H_tb, r_tb					)
		call kspace_allocator(		H_tb, r_tb, 			en_k, V_ka, A_ka, Om_kab	)
		call print_basis_info()
		!
		!----------------------------------------------------------------------------------------------------------------------------------
		!	loop	k-space
		!----------------------------------------------------------------------------------------------------------------------------------
		do kiz = 1, mp_grid(3)
			do kiy = 1, mp_grid(2)
				do kix = 1, mp_grid(1)
					!
					!
					ki	=	get_rel_kpt(kix, kiy, kiz, kpt)
					if( mpi_ki_selector(ki, num_kpts)	) then
						!----------------------------------------------------------------------------------------------------------------------------------
						!----------------------------------------------------------------------------------------------------------------------------------
						!----------------------------------------------------------------------------------------------------------------------------------
						!	INTERPOLATE
						!----------------------------------------------------------------------------------------------------------------------------------
						call get_wann_interp(do_gauge_trafo, H_tb, r_tb, a_latt, recip_latt, R_vect, ki, kpt(:), 	en_k, V_ka, A_ka, Om_kab )
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	MEP
						!----------------------------------------------------------------------------------------------------------------------------------
						if(	allocated(mep_tens_ic_loc)	)		mep_tens_ic_loc	=	mep_tens_ic_loc + 	mep_niu_IC(V_ka, en_k)		!	itinerant		(Kubo)
						if( allocated(mep_tens_lc_loc)	)		mep_tens_lc_loc	=	mep_tens_lc_loc + 	mep_niu_LC(V_ka, en_k)		!	local			(Kubo)
						if( allocated(mep_tens_cs_loc)	)		mep_tens_cs_loc =	mep_tens_cs_loc + 	mep_niu_CS(A_ka, Om_kab)	!	chern simons	(geometrical)
						!
						if( allocated(mep_2014_loc)		)		mep_2014_loc	=	mep_2014_loc	+	mep_niu2014_full(V_ka, en_k)
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	KUBO MEP (MEP with fermi_dirac)
						!----------------------------------------------------------------------------------------------------------------------------------
						if( allocated(kubo_mep_ic_loc)	)		kubo_mep_ic_loc	=	kubo_mep_ic_loc +	kubo_mep_IC(eFermi, T_kelvin, V_ka, en_k, ic_skipped)
						if( allocated(kubo_mep_lc_loc)	)		kubo_mep_lc_loc	=	kubo_mep_lc_loc +	kubo_mep_LC(eFermi, T_kelvin, V_ka, en_k, lc_skipped)
						if( allocated(kubo_mep_cs_loc)	)		kubo_mep_cs_loc	=	kubo_mep_cs_loc +	kubo_mep_CS(eFermi, T_kelvin, en_k, A_ka, Om_kab)
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	AHC
						!----------------------------------------------------------------------------------------------------------------------------------
						if(allocated(kubo_ahc_loc)		)		kubo_ahc_loc	=	kubo_ahc_loc	+ 	kubo_ahc_tens(en_k,	Om_kab, eFermi, T_kelvin, unit_vol)
						if(allocated(velo_ahc_loc)		)		velo_ahc_loc	=	velo_ahc_loc	+	velo_ahc_tens(en_k, V_ka,	eFermi, T_kelvin, unit_vol)
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	OPT
						!----------------------------------------------------------------------------------------------------------------------------------
						if(allocated(kubo_opt_s_loc) .and. allocated(kubo_opt_a_loc)	)		then
							call kubo_opt_tens(hw, unit_vol, eFermi, T_kelvin, i_eta_smr, en_k, A_ka, 		tempS, tempA)
							kubo_opt_s_loc	=	kubo_opt_s_loc	+	tempS							
							kubo_opt_a_loc	=	kubo_opt_a_loc	+	tempA
						end if
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	GYRO
						!----------------------------------------------------------------------------------------------------------------------------------
						if(	allocated(gyro_C_loc)		)		gyro_C_loc		=	gyro_C_loc		+	get_gyro_C(en_k, V_ka, eFermi, T_kelvin)
						if(	allocated(gyro_D_loc)		)		gyro_D_loc		=	gyro_D_loc		+ 	get_gyro_D(en_k, V_ka, Om_kab, eFermi, T_kelvin)
						if(	allocated(gyro_Dw_loc)		)		gyro_Dw_loc		=	gyro_Dw_loc		+ 	get_gyro_Dw()	!dummy returns 0 currently!!!
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	$WILDCARD	:	add new tensors here!
						!----------------------------------------------------------------------------------------------------------------------------------

						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	COUNT ELECTRONS
						!----------------------------------------------------------------------------------------------------------------------------------
						call fd_count_el(en_k, eFermi, T_kelvin,		N_el_k, sum_N_el_loc, min_n_el, max_n_el)
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!----------------------------------------------------------------------------------------------------------------------------------
						!----------------------------------------------------------------------------------------------------------------------------------
						!
						n_ki_loc = n_ki_loc + 1
					end if
					!
					!
				end do
			end do
		end do
		!
		!
		!----------------------------------------------------------------------------------------------------------------------------------
		!	PRINT INFO
		!----------------------------------------------------------------------------------------------------------------------------------
		call print_core_info(n_ki_loc, sum_N_el_loc, min_n_el, max_n_el)
		!
		!----------------------------------------------------------------------------------------------------------------------------------
		!	REDUCE MPI & WRITE FILES
		!----------------------------------------------------------------------------------------------------------------------------------
		call norm_K_int_and_write(			n_ki_loc,																	&
										mep_2014_loc,																	&
										mep_tens_ic_loc, mep_tens_lc_loc, mep_tens_cs_loc,								&	
										kubo_mep_ic_loc, kubo_mep_lc_loc, kubo_mep_cs_loc,								&									
										kubo_ahc_loc, velo_ahc_loc,														&
										kubo_opt_s_loc, kubo_opt_a_loc,													&
										gyro_C_loc, gyro_D_loc, gyro_Dw_loc												&
							)
		!----------------------------------------------------------------------------------------------------------------------------------
		!----------------------------------------------------------------------------------------------------------------------------------
		!----------------------------------------------------------------------------------------------------------------------------------	
		return
	end subroutine






!private:
	!
	!
!
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!	HELPERS
	!----------------------------------------------------------------------------------------------------------------------------------	
	subroutine norm_K_int_and_write(		n_ki_loc, 																&
									mep_2014_loc,																		&
									mep_tens_ic_loc, mep_tens_lc_loc, mep_tens_cs_loc,									&	
									kubo_mep_ic_loc, kubo_mep_lc_loc, kubo_mep_cs_loc,									&									
									kubo_ahc_loc, velo_ahc_loc,															&
									kubo_opt_s_loc, kubo_opt_a_loc,														&
									gyro_C_loc, gyro_D_loc, gyro_Dw_loc													&
							)
		!-----------------------------------------------------------------------------------------------------------
		!		local:		sum of k-points	 over local( within single mpi) thread
		!		global:		sum of k-points over whole mesh
		!
		!			mep_tens_ic_loc, , mep_tens_lc_loc, mep_tens_cs_loc ar given as arrays over valence bands
		!-----------------------------------------------------------------------------------------------------------
		integer,						intent(in)		::	n_ki_loc 	
		real(dp),		allocatable,	intent(inout)	::	mep_2014_loc(:,:),														&
															mep_tens_ic_loc(:,:,:), mep_tens_lc_loc(:,:,:), mep_tens_cs_loc(:,:,:),	&
															kubo_mep_ic_loc(:,:), kubo_mep_lc_loc(:,:), kubo_mep_cs_loc(:,:),		&				
															kubo_ahc_loc(:,:), velo_ahc_loc(:,:)									
		complex(dp),	allocatable,	intent(inout)	::	kubo_opt_s_loc(:,:), kubo_opt_a_loc(:,:),								&
															gyro_C_loc(:,:), gyro_D_loc(:,:), gyro_Dw_loc(:,:)		
		!-----------------------------------------------------------------------------------------------------------
		integer											::	n_ki_glob
		real(dp), 				allocatable				::	mep_bands_loc(:,:,:),													&
															mep_bands_glob(:,:,:),													&
															mep_2014_glob(:,:),														&
															!
															mep_sum_ic_loc(:,:), mep_sum_lc_loc(:,:), mep_sum_cs_loc(:,:),			&
															mep_sum_ic_glob(:,:), mep_sum_lc_glob(:,:), mep_sum_cs_glob(:,:),		&
															kubo_mep_ic_glob(:,:), kubo_mep_lc_glob(:,:), kubo_mep_cs_glob(:,:),	&
															!
															kubo_ahc_glob(:,:), velo_ahc_glob(:,:)
		complex(dp),			allocatable				::	kubo_opt_s_glob(:,:), kubo_opt_a_glob(:,:),								&
															gyro_C_glob(:,:), gyro_D_glob(:,:), gyro_Dw_glob(:,:)		
		!
		!-----------------------------------------------------------------------------------------------------------
		!			sum over bands before communication	(else: expensive )
		!-----------------------------------------------------------------------------------------------------------		!
		call mep_sum_bands(		mep_tens_ic_loc, 	mep_tens_lc_loc, 	mep_tens_cs_Loc,	&
								mep_sum_ic_loc,		mep_sum_lc_loc,		mep_sum_cs_loc,		&
								mep_bands_loc												&
			)
		!
		!-----------------------------------------------------------------------------------------------------------
		!			COMMUNICATION (k-pts)
		!-----------------------------------------------------------------------------------------------------------
		write(*,'(a,i3,a)') "[#",mpi_id,"; core_worker]:  now start reduction"
		!
		!	***********				SUM OVER NODES				******************************************************
		call mpi_reduce_tens(	n_ki_loc,				n_ki_glob		)
		!
		call mpi_reduce_tens(	mep_bands_loc	,	mep_bands_glob		)
		call mpi_reduce_tens(	mep_2014_loc	,	mep_2014_glob		)		
		call mpi_reduce_tens(	mep_sum_ic_loc	,	mep_sum_ic_glob		)
		call mpi_reduce_tens(	mep_sum_lc_loc	,	mep_sum_lc_glob		)
		call mpi_reduce_tens(	mep_sum_cs_loc	,	mep_sum_cs_glob		)
		call mpi_reduce_tens(	kubo_mep_ic_loc	,	kubo_mep_ic_glob	)
		call mpi_reduce_tens(	kubo_mep_lc_loc	,	kubo_mep_lc_glob	)
		call mpi_reduce_tens(	kubo_mep_cs_loc	,	kubo_mep_cs_glob	)
		call mpi_reduce_tens(	kubo_ahc_loc	,	kubo_ahc_glob		)
		call mpi_reduce_tens(	velo_ahc_loc	,	velo_ahc_glob		)
		call mpi_reduce_tens(	kubo_opt_s_loc	,	kubo_opt_s_glob		)
		call mpi_reduce_tens(	kubo_opt_a_loc	,	kubo_opt_a_glob		)
		call mpi_reduce_tens(	gyro_C_loc		,	gyro_C_glob			)
		call mpi_reduce_tens(	gyro_D_loc		,	gyro_D_glob			)
		call mpi_reduce_tens(	gyro_Dw_loc		,	gyro_Dw_glob		)
		!
		if(mpi_id == mpi_root_id)	then
			write(*,*)				""
			write(*,*)				""
			write(*,*)				"..."					
			write(*,'(a,i3,a,i8,a)') "[#",mpi_id,"; core_worker]: collected tensors from",mpi_nProcs," mpi-threads"
		end if
		!-----------------------------------------------------------------------------------------------------------
		!			NORMALIZE K-SPACE INTEGRAL
		!-----------------------------------------------------------------------------------------------------------
		call normalize_k_int(mep_2014_glob)
		call normalize_k_int(mep_sum_ic_glob)
		call normalize_k_int(mep_sum_lc_glob)
		call normalize_k_int(mep_sum_cs_glob)
		call normalize_k_int(mep_bands_glob)
		call normalize_k_int(kubo_ahc_glob)
		call normalize_k_int(velo_ahc_glob)
		call normalize_k_int(kubo_opt_s_glob)
		call normalize_k_int(kubo_opt_a_glob)
		call normalize_k_int(gyro_C_glob)
		call normalize_k_int(gyro_D_glob)
		call normalize_k_int(gyro_Dw_glob)
		!
		!-----------------------------------------------------------------------------------------------------------
		!			OUTPUT
		!-----------------------------------------------------------------------------------------------------------
		if(mpi_id == mpi_root_id) then
			write(*,*)	"*"
			write(*,*)	"------------------OUTPUT----------------------------------------------"
			call write_mep_bands(			n_ki_glob,		mep_bands_glob												)
			call write_mep_tensors(			n_ki_glob,		mep_2014_glob ,											&
															mep_sum_ic_glob, 	mep_sum_lc_glob, 	mep_sum_cs_glob		)
			call write_kubo_mep_tensors(	n_ki_glob,		kubo_mep_ic_glob, 	kubo_mep_lc_glob, 	kubo_mep_cs_glob	)
			call write_ahc_tensor(			n_ki_glob,		kubo_ahc_glob, 		velo_ahc_glob							)
			call write_opt_tensors(			n_ki_glob,		kubo_opt_s_glob, 	kubo_opt_a_glob							)			
			call write_gyro_tensors( 		n_ki_glob,		gyro_C_glob, 		gyro_D_glob, 		gyro_Dw_glob		)
			write(*,*)	"*"
			write(*,*)	"----------------------------------------------------------------"
		end if	
		!
		!
		!-----------------------------------------------------------------------------------------------------------
		!			DEBUG
		!-----------------------------------------------------------------------------------------------------------
		if (mpi_id	==	mpi_root_id)  call print_debug_info(n_ki_glob)
		!
		!
		return
	end subroutine








	subroutine mep_sum_bands(		mep_tens_ic_loc, 	mep_tens_lc_loc, 	mep_tens_cs_Loc,	&
									mep_sum_ic_loc,		mep_sum_lc_loc,		mep_sum_cs_loc,		&
									mep_bands_loc						&
								)
		real(dp),		allocatable, 		intent(inout)		::	mep_tens_ic_loc(:,:,:), mep_tens_lc_loc(:,:,:), mep_tens_cs_loc(:,:,:),		& 
																	mep_sum_ic_loc(:,:), 	mep_sum_lc_loc(:,:), 	mep_sum_cs_loc(:,:),		& 
																	mep_bands_loc(:,:,:)
		integer													::	n0
		!
		if(	allocated(mep_tens_ic_loc) .and. allocated(mep_tens_lc_loc) .and. allocated(mep_tens_cs_Loc)	) then
			if(.not. allocated(	mep_sum_ic_loc	))	allocate(	mep_sum_ic_Loc(		3,	3					))		
			if(.not. allocated(	mep_sum_lc_loc	))	allocate(	mep_sum_lc_Loc(		3,	3					))
			if(.not. allocated(	mep_sum_cs_loc	))	allocate(	mep_sum_cs_Loc(		3,	3					))
			if(.not. allocated(	mep_bands_loc	))	allocate(	mep_bands_loc(		3,	3,	valence_bands	))
			!
			mep_sum_ic_loc		=	0.0_dp
			mep_sum_lc_loc		=	0.0_dp
			mep_sum_cs_loc		=	0.0_dp
			mep_bands_loc		=	0.0_dp	
			do n0 = 1, valence_bands
				mep_sum_ic_loc(:,:)	=	mep_sum_ic_loc(:,:)	+	mep_tens_ic_loc(:,:,n0)
				mep_sum_lc_loc(:,:)	=	mep_sum_lc_loc(:,:)	+	mep_tens_lc_loc(:,:,n0)
				mep_sum_cs_loc(:,:)	=	mep_sum_cs_loc(:,:)	+	mep_tens_cs_loc(:,:,n0)
				!
				if(do_write_mep_bands) then
					mep_bands_loc(:,:,n0)	=	mep_tens_ic_loc(:,:,n0) + mep_tens_lc_loc(:,:,n0)	+ mep_tens_cs_loc(:,:,n0)
				end if			
			end do
			!
		end if
		!
		return 
	end subroutine
!
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!	ALLOCATORS
	!----------------------------------------------------------------------------------------------------------------------------------
	subroutine allo_core_loc_arrays(			N_el_k, sum_N_el_loc, max_n_el, min_n_el,										&
										mep_2014_loc,																	&
										mep_tens_ic_loc, mep_tens_lc_loc, mep_tens_cs_loc,								&	
										kubo_mep_ic_loc, kubo_mep_lc_loc, kubo_mep_cs_loc,								&									
										kubo_ahc_loc, velo_ahc_loc,														&
										tempS, tempA, kubo_opt_s_loc, kubo_opt_a_loc,									&
										gyro_C_loc, gyro_D_loc, gyro_Dw_loc												&
							)
		real(dp),						intent(inout)	::	N_el_k, sum_N_el_loc, max_n_el, min_n_el
		real(dp),		allocatable,	intent(inout)	::	mep_2014_loc(:,:),														&
															mep_tens_ic_loc(:,:,:), mep_tens_lc_loc(:,:,:), mep_tens_cs_loc(:,:,:),	&
															kubo_mep_ic_loc(:,:), kubo_mep_lc_loc(:,:), kubo_mep_cs_loc(:,:),		&				
															kubo_ahc_loc(:,:), velo_ahc_loc(:,:)									
		complex(dp),	allocatable,	intent(inout)	::	tempS(:,:), tempA(:,:), kubo_opt_s_loc(:,:), kubo_opt_a_loc(:,:),		&
															gyro_C_loc(:,:), gyro_D_loc(:,:), gyro_Dw_loc(:,:)		
		!----------------------------------------------------------------------------------------------------------------------------------
		!	ALLOCATE & INIT
		!----------------------------------------------------------------------------------------------------------------------------------													
		N_el_k				=	0.0_dp
		sum_N_el_loc		=	0.0_dp
		max_n_el			=	0.0_dp
		min_n_el			=	2.0_dp
		!
		if( do_mep	)	then
			allocate(	mep_tens_ic_loc(	3,3,	valence_bands	)	)
			allocate(	mep_tens_lc_loc(	3,3,	valence_bands	)	)
			allocate(	mep_tens_cs_loc(	3,3,	valence_bands	)	)
			allocate(	mep_2014_loc(				3,3				)	)
			!
			mep_2014_loc		=	0.0_dp
			mep_tens_ic_loc		=	0.0_dp
			mep_tens_lc_loc		=	0.0_dp
			mep_tens_cs_loc		=	0.0_dp
		end if
		!
		if(	do_ahc	)	then		
			allocate(	kubo_ahc_loc(				3,3	)			)	
			allocate(	velo_ahc_loc(				3,3	)			)
			!
			kubo_ahc_loc		=	0.0_dp
			velo_ahc_loc		=	0.0_dp	
		end if
		!
		if(	do_kubo	)	then	
			allocate(	kubo_mep_ic_loc(			3,3	)			)	
			allocate(	kubo_mep_lc_loc(			3,3	)			)	
			allocate(	kubo_mep_cs_loc(			3,3	)			)
			!
			kubo_mep_ic_loc		=	0.0_dp
			kubo_mep_lc_loc		=	0.0_dp
			kubo_mep_cs_loc		=	0.0_dp
		end if
		!
		if(	do_opt	)	then
			allocate(	tempS(						3,3	)			)
			allocate(	tempA(						3,3	)			)				
			allocate(	kubo_opt_s_loc(				3,3	)			)
			allocate(	kubo_opt_a_loc(				3,3	)			)
			!
			kubo_opt_a_loc		=	0.0_dp
			kubo_opt_s_loc		=	0.0_dp
		end if
		!
		if(	do_gyro	)	then
			allocate(	gyro_D_loc(					3,3	)			)
			allocate(	gyro_Dw_loc(				3,3	)			)																								
			allocate(	gyro_C_loc(					3,3	)			)
			!
			gyro_D_loc			=	0.0_dp
			gyro_Dw_loc			=	0.0_dp
			gyro_C_loc			=	0.0_dp	
		end if	
		!
		!
		return
	end subroutine



!
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!	PRINTERS
	!----------------------------------------------------------------------------------------------------------------------------------
	subroutine print_basis_info()
		if(use_mpi)	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		if(mpi_id == mpi_root_id) then
			write(*,*)	"*"
			write(*,*)	"----------------------------------------------------------------"
			write(*,*)	"----------------------------------------------------------------"
			write(*,*)	"----------------------------------------------------------------"
			write(*,*)	"*"
			write(*,*)	"*"
			if(do_gauge_trafo)	then
					write(*,*)	"***^^^^	-	FULL INTERPOLATION MODE	-	^^^^***"
			else
					write(*,*)	"***^^^^	-	WANNIER GAUGE MODE 	-	^^^^***"
			end if
			write(*,*)	"*"
			write(*,*)	"*"
			write(*,*)	"----------------------------------------------------------------"
			write(*,*)	"----------------------------------------------------------------"
			write(*,*)	"----------------------------------------------------------------"
		end if
		if(use_mpi)	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		write(*,'(a,i3,a,a,a,i4,a)')		"[#",mpi_id,"; core_worker/",cTIME(time()),		&
											"]:  I start interpolating now (nValence=",valence_bands,")...."
		!
		!
		return
	end subroutine
	!
	!
	subroutine print_core_info(n_ki_loc, sum_N_el_loc, min_n_el, max_n_el)
		integer,		intent(in)		::	n_ki_loc
		real(dp),		intent(in)		::	sum_N_el_loc, min_n_el, max_n_el
		character(len=60)				::	gauge_label
		!
		if(do_gauge_trafo)	then
			write(gauge_label,*)		'in the (H) Hamiltonian gauge'
		else
			write(gauge_label,*)		'in the (W) Wannier gauge (no gauge trafo performed!)'
		end if
		!
		!
		write(*,'(a,i3,a,a,a,i8,a)')					"[#",mpi_id,"; core_worker/",cTIME(time()),		&
													"]: ...finished interpolating ",n_ki_loc," kpts "//gauge_label
		!
		if(use_mpi)		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		!
		!
		if(n_ki_loc > 0)	then
			write(*,'(a,i3,a,f5.2,a,f5.2,a,f5.2,a)')	"[#",mpi_id,"; core_worker]: avg el count ",						&
																						sum_N_el_loc/real(n_ki_loc,dp),		&
																					"	(min: ", 							&
																						min_n_el,							&
																					";   max: ", 							&
																						max_n_el,							&
																					")"
		end if
		!
		!
		return
	end subroutine
	!
	!
	subroutine print_debug_info(n_ki_glob)
		integer,		intent(in)		::	n_ki_glob
		!
		if(	n_ki_glob	/=	num_kpts .or. n_ki_glob <= 0	) then
				write(*,'(a,i3,a,i12)')		"[#",mpi_id,";core_worker]	ERROR n_ki_glob=",n_ki_glob
				write(*,'(a,i3,a,i12)')		"[#",mpi_id,";core_worker]	ERROR  num_kpts=",num_kpts
			stop "[core_worker]:	ERROR n_ki_glob is not equal to the given mp_grid"
		end if
		!
		!ToDo:	compare  mep_sum vs sum {mep_bands}

		return
	end subroutine


	!
	!
end module core


