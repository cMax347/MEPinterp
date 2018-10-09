module core
	!this module uses a semiclassic approach to calculate the first ordrer correction
	!	to the polariztion induced by a perturbive magnetic field
	! 	see Niu PRL 112, 166601 (2014)
	!use omp_lib
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
	use mpi
	use matrix_math,	only:	my_Levi_Civita
	use constants,		only:	dp, aUtoAngstrm, auToTesla,						&			
								mpi_root_id, mpi_id, mpi_nProcs, ierr			
	use statistics,		only:	fd_get_N_el
	use input_paras,	only:	a_latt, 										&
								valence_bands, 									&
								seed_name,										&
								kubo_tol,										&
								hw, eFermi, T_kelvin, i_eta_smr, 				&
								unit_vol
	!
	use k_space,		only:	get_recip_latt, get_mp_grid, 					&
								get_rel_kpt,									&
								normalize_k_int
	use wann_interp,	only:	get_wann_interp
	!
	use mep_niu,		only:	mep_niu_CS,	mep_niu_IC, mep_niu_LC
	use kubo,			only:	kubo_ahc_tens, velo_ahc_tens, kubo_opt_tens
	use kubo_mep,		only:	kubo_mep_CS, kubo_mep_LC, kubo_mep_IC
	use gyro,			only:	get_gyro_C, get_gyro_D, get_gyro_Dw
	!
	use file_io,		only:	mpi_read_tb_basis,								&
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
												kpt(3),	recip_latt(3,3),				&
												!local sum targets:
												mep_tens_ic_loc(	3,3),	 			&
												mep_tens_lc_loc(	3,3),				&
												mep_tens_cs_loc(	3,3),				&
												kubo_ahc_loc(		3,3),				&
												velo_ahc_loc(		3,3),				&
												kubo_mep_ic_loc(	3,3),				&
												kubo_mep_lc_loc(	3,3),				&
												kubo_mep_cs_loc(	3,3)
												!
		complex(dp)							::	tempS(3,3), 	tempA(3,3),				&
												kubo_opt_s_loc(		3,3),				&
												kubo_opt_a_loc(		3,3),				&
												gyro_D_loc(			3,3),				&
												gyro_Dw_loc(		3,3),				&																								
												gyro_C_loc(			3,3)																								
												!
		integer								::	kix, kiy, kiz, ki, 						&
												mp_grid(3), n_ki_loc, 					&
												ic_skipped, lc_skipped
		complex(dp),	allocatable			::	H_tb(:,:,:), r_tb(:,:,:,:), 			&
												A_ka(:,:,:), Om_kab(:,:,:,:),			&
												V_ka(:,:,:)
		real(dp),		allocatable			::	en_k(:), R_vect(:,:)
		!
		!----------------------------------------------------------------------------------------------------------------------------------
		!	INIT
		!----------------------------------------------------------------------------------------------------------------------------------
		N_el_k				=	0.0_dp
		sum_N_el_loc		=	0.0_dp
		max_n_el			=	0.0_dp
		min_n_el			=	2.0_dp
		mep_tens_ic_loc		=	0.0_dp
		mep_tens_lc_loc		=	0.0_dp
		mep_tens_cs_loc		=	0.0_dp
		kubo_mep_ic_loc		=	0.0_dp
		kubo_mep_lc_loc		=	0.0_dp
		kubo_mep_cs_loc		=	0.0_dp
		kubo_ahc_loc		=	0.0_dp
		kubo_opt_a_loc		=	0.0_dp
		kubo_opt_s_loc		=	0.0_dp
		gyro_D_loc			=	0.0_dp
		gyro_Dw_loc			=	0.0_dp
		gyro_C_loc			=	0.0_dp
		!
		!
		!----------------------------------------------------------------------------------------------------------------------------------
		!	get 	k-space
		!----------------------------------------------------------------------------------------------------------------------------------
		mp_grid				=	get_mp_grid()
		n_ki_loc			= 	0
		num_kpts			= 	mp_grid(1)*mp_grid(2)*mp_grid(3)
		recip_latt			=	get_recip_latt()
		write(*,*)''
		!
		!----------------------------------------------------------------------------------------------------------------------------------
		!	get		TB	BASIS
		!----------------------------------------------------------------------------------------------------------------------------------
		call mpi_read_tb_basis(	seed_name, R_vect,		 	H_tb, r_tb					)
		call kspace_allocator(		H_tb, r_tb, 			en_k, V_ka, A_ka, Om_kab	)
		!
		!
		!
		write(*,*)''
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		write(*,'(a,i3,a,a,a,i4,a)')		"[#",mpi_id,"; core_worker/",cTIME(time()),		&
											"]:  I start interpolating now (nValence=",valence_bands,")...."
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
						call get_wann_interp(H_tb, r_tb, a_latt, recip_latt, R_vect, kpt(:), 	en_k, V_ka, A_ka, Om_kab )
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	MEP
						!----------------------------------------------------------------------------------------------------------------------------------
						mep_tens_ic_loc	=	mep_tens_ic_loc + 	mep_niu_IC(V_ka, en_k)		!	itinerant		(Kubo)
						mep_tens_lc_loc	=	mep_tens_lc_loc + 	mep_niu_LC(V_ka, en_k)		!	local			(Kubo)
						mep_tens_cs_loc =	mep_tens_cs_loc + 	mep_niu_CS(A_ka, Om_kab)	!	chern simons	(geometrical)
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	KUBO MEP (MEP with fermi_dirac)
						!----------------------------------------------------------------------------------------------------------------------------------
						kubo_mep_ic_loc	=	kubo_mep_ic_loc +	kubo_mep_IC(eFermi, T_kelvin, V_ka, en_k, ic_skipped)
						kubo_mep_lc_loc	=	kubo_mep_lc_loc +	kubo_mep_LC(eFermi, T_kelvin, V_ka, en_k, lc_skipped)
						kubo_mep_cs_loc	=	kubo_mep_cs_loc +	kubo_mep_CS(eFermi, T_kelvin, en_k, A_ka, Om_kab)
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	AHC
						!----------------------------------------------------------------------------------------------------------------------------------
						kubo_ahc_loc	=	kubo_ahc_loc	+ 	kubo_ahc_tens(en_k,	Om_kab, eFermi, T_kelvin, unit_vol)
						velo_ahc_loc	=	velo_ahc_loc	+	velo_ahc_tens(en_k, V_ka,	eFermi, T_kelvin, unit_vol)
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	OPT
						!----------------------------------------------------------------------------------------------------------------------------------
						call kubo_opt_tens(hw, unit_vol, eFermi, T_kelvin, i_eta_smr, en_k, A_ka, 		tempS, tempA)
						kubo_opt_s_loc	=	kubo_opt_s_loc	+	tempS							
						kubo_opt_a_loc	=	kubo_opt_a_loc	+	tempA
						!
						!----------------------------------------------------------------------------------------------------------------------------------
						!	GYRO
						!----------------------------------------------------------------------------------------------------------------------------------
						gyro_C_loc		=	gyro_C_loc		+	get_gyro_C(en_k, V_ka, eFermi, T_kelvin)
						gyro_D_loc		=	gyro_D_loc		+ 	get_gyro_D(en_k, V_ka, Om_kab, eFermi, T_kelvin)
						gyro_Dw_loc		=	gyro_Dw_loc		+ 	get_gyro_Dw()	!dummy returns 0 currently!!!
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
		call sumK_and_print(			n_ki_loc, sum_N_el_loc,															&
										mep_tens_ic_loc, mep_tens_lc_loc, mep_tens_cs_loc,								&	
										kubo_mep_ic_loc, kubo_mep_lc_loc, kubo_mep_cs_loc,								&									
										kubo_ahc_loc, velo_ahc_loc,														&
										kubo_opt_s_loc, kubo_opt_a_loc,													&
										gyro_C_loc, gyro_D_loc, gyro_Dw_loc												&
							)
		!
		!----------------------------------------------------------------------------------------------------------------------------------
		!----------------------------------------------------------------------------------------------------------------------------------
		!----------------------------------------------------------------------------------------------------------------------------------	
		return
	end subroutine








!private:
	!
	!
!HELPERS
	subroutine sumK_and_print(		n_ki_loc, sum_N_el_loc,																&
									mep_tens_ic_loc, mep_tens_lc_loc, mep_tens_cs_loc,									&	
									kubo_mep_ic_loc, kubo_mep_lc_loc, kubo_mep_cs_loc,									&									
									kubo_ahc_loc, velo_ahc_loc,															&
									kubo_opt_s_loc, kubo_opt_a_loc,														&
									gyro_C_loc, gyro_D_loc, gyro_Dw_loc													&
							)
		!		local:		sum of k-points	 over local( within single mpi) thread
		!		global:		sum of k-points over whole mesh
		!
		integer,			intent(in)		::	n_ki_loc 	
		real(dp),			intent(in)		::	sum_N_el_loc,															&
												mep_tens_ic_loc(3,3), mep_tens_lc_loc(3,3), mep_tens_cs_loc(3,3),		&
												kubo_mep_ic_loc(3,3), kubo_mep_lc_loc(3,3), kubo_mep_cs_loc(3,3),		&
												kubo_ahc_loc(3,3), velo_ahc_loc(3,3)														
		complex(dp),		intent(in)		::	kubo_opt_s_loc(3,3), kubo_opt_a_loc(3,3),								&
												gyro_C_loc(3,3), gyro_D_loc(3,3), gyro_Dw_loc(3,3)		
		!								---------------------------------------------------------------
		integer								::	n_ki_glob
		!								---------------------------------------------------------------
		real(dp)							::	avg_N_el_glob,															&
												mep_tens_ic_glob(3,3), mep_tens_lc_glob(3,3), mep_tens_cs_glob(3,3),	&
												kubo_mep_ic_glob(3,3), kubo_mep_lc_glob(3,3), kubo_mep_cs_glob(3,3),	&
												kubo_ahc_glob(3,3), velo_ahc_glob(3,3)
		complex(dp)							::	kubo_opt_s_glob(3,3), kubo_opt_a_glob(3,3),								&
												gyro_C_glob(3,3), gyro_D_glob(3,3), gyro_Dw_glob(3,3)		
		!								---------------------------------------------------------------
		avg_N_el_glob		=	0.0_dp								
		!	
		mep_tens_ic_glob	=	0.0_dp
		mep_tens_lc_glob	=	0.0_dp
		mep_tens_cs_glob	=	0.0_dp	
		kubo_mep_ic_glob	=	0.0_dp
		kubo_mep_lc_glob	=	0.0_dp
		kubo_mep_cs_glob	=	0.0_dp	
		kubo_ahc_glob		=	0.0_dp
		velo_ahc_glob		=	0.0_dp
		kubo_opt_a_glob		=	0.0_dp
		kubo_opt_s_glob		=	0.0_dp
		gyro_C_glob			=	0.0_dp
		gyro_D_glob			=	0.0_dp
		gyro_Dw_glob		=	0.0_dp
		!
		!
		!sum mep over global kpts
		call MPI_REDUCE(	n_ki_loc,				n_ki_glob,			1,		MPI_INTEGER,		MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	sum_N_el_loc,		avg_N_el_glob,			1,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		!
		call MPI_REDUCE(	mep_tens_ic_loc,	mep_tens_ic_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		, 	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	mep_tens_lc_loc,	mep_tens_lc_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		, 	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	mep_tens_cs_loc,	mep_tens_cs_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		!
		call MPI_REDUCE(	kubo_mep_ic_loc,	kubo_mep_ic_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	kubo_mep_lc_loc,	kubo_mep_lc_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	kubo_mep_cs_loc,	kubo_mep_cs_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		!
		call MPI_REDUCE(	kubo_ahc_loc,		kubo_ahc_glob,			9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	velo_ahc_loc,		velo_ahc_glob,			9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)		
		!
		call MPI_REDUCE(	kubo_opt_s_loc,		kubo_opt_s_glob,		9,	MPI_DOUBLE_COMPLEX,		MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD,	ierr)
		call MPI_REDUCE(	kubo_opt_a_loc,		kubo_opt_a_glob,		9,	MPI_DOUBLE_COMPLEX,		MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD,	ierr)
		!
		call MPI_REDUCE(	gyro_C_loc,			gyro_C_glob,			9,	MPI_DOUBLE_COMPLEX,		MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	gyro_D_loc,			gyro_D_glob,			9,	MPI_DOUBLE_COMPLEX,		MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	gyro_Dw_loc,		gyro_Dw_glob,			9,	MPI_DOUBLE_COMPLEX,		MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		!
		!
		!write tensors to files
		if(mpi_id == mpi_root_id) 	then
			write(*,*)	"--------------------------------------------------------------------------------------------------------"	
			if(	n_ki_glob	/=	num_kpts	) stop "[core_worker]:	ERROR n_ki_glob is not equal to the given mp_grid"
			!
			!
			!	NORMALIZE INTEGRATION
			avg_N_el_glob	=	avg_N_el_glob	/	real(n_ki_glob,dp)
			call normalize_k_int(mep_tens_ic_glob)
			call normalize_k_int(mep_tens_lc_glob)
			call normalize_k_int(mep_tens_cs_glob)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; core_worker]: calculated MEP tensor on ",n_ki_glob," kpts"
			!
			call normalize_k_int(kubo_mep_ic_glob)
			call normalize_k_int(kubo_mep_lc_glob)
			call normalize_k_int(kubo_mep_cs_glob)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; core_worker]: calculated  KUBO MEP tensor on ",n_ki_glob," kpts"
			!
			!
			kubo_ahc_glob	=	kubo_ahc_glob		/	real(n_ki_glob,dp)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; core_worker]: calculated AHC tensor on ",n_ki_glob," kpts"
			!			!
			velo_ahc_glob	=	velo_ahc_glob		/	real(n_ki_glob,dp)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; core_worker]: calculated AHC (via velo) tensor on ",n_ki_glob," kpts"
			!
			kubo_opt_s_glob	=	kubo_opt_s_glob		/	real(n_ki_glob,dp)			
			kubo_opt_a_glob	=	kubo_opt_a_glob		/	real(n_ki_glob,dp)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; core_worker]: calculated OPT tensor on ",n_ki_glob," kpts"
			!
			call normalize_k_int(gyro_C_glob)
			call normalize_k_int(gyro_D_glob)			
			call normalize_k_int(gyro_Dw_glob)			

			!
			!	WRITE FILES
			call write_mep_tensors(mep_tens_ic_glob, mep_tens_lc_glob, mep_tens_cs_glob)
			call write_kubo_mep_tensors(kubo_mep_ic_glob, kubo_mep_lc_glob, kubo_mep_cs_glob)
			call write_ahc_tensor(kubo_ahc_glob, velo_ahc_glob)
			call write_opt_tensors(kubo_opt_s_glob, kubo_opt_a_glob)
			call write_gyro_tensors( gyro_C_glob, gyro_D_glob, gyro_Dw_glob)
			write(*,*)	""
			write(*,*)	"--------------------------------------------------------------------------------------------------------"	
		end if
		!
		return
	end subroutine





	pure subroutine fd_count_el(en_k, eFermi, T_kelvin, el_count, sum_loc, n_el_min, n_el_max)
		real(dp),		intent(in)					::		en_k(:), eFermi, T_kelvin
		real(dp),		intent(out)					::		el_count
		real(dp),		intent(inout)				::		sum_loc, n_el_min, n_el_max
		!
		el_count			=	fd_get_N_el(en_k, eFermi, T_kelvin)
		sum_loc				=	sum_loc		+	el_count
		!		!
		if(	el_count	<	n_el_min	)		n_el_min	= 	el_count
		if( el_count	>	n_el_max	)		n_el_max	=	el_count
		!
		return
	end subroutine


	subroutine kspace_allocator(H_tb, r_tb, en_k, V_ka, A_ka, Om_kab)
		complex(dp),	allocatable, intent(inout)	::		H_tb(:,:,:), r_tb(:,:,:,:),					& 
															V_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp),		allocatable, intent(inout)	::		en_k(:)
		!
		!allocate k-space
		allocate(	en_k(						size(H_tb,2)	)	)
		allocate(	V_ka(	3,	size(H_tb,1),	size(H_tb,2)	)	)
		!
		if(	allocated(r_tb)	)	then	
			allocate(	A_ka(	  3,		size(r_tb,2),	size(r_tb,3)	)	)
			allocate(	Om_kab(	3,	3,		size(r_tb,2),	size(r_tb,3)	)	)
			write(*,'(a,i3,a)')		"[#",mpi_id,"; core_worker]: will use position operator"	
		end if
		!
		return
	end subroutine


	logical pure function mpi_ki_selector(ki_request, num_kpts)
		integer,		intent(in)		::		ki_request, num_kpts 
		integer							::		ki_todo
		!
		mpi_ki_selector = .false.
		!
		if(		mpi_nProcs == 1		) then
			mpi_ki_selector	=	.true.
		else
			loop_todos: do ki_todo = mpi_id +1, num_kpts, mpi_nProcs
				if(	ki_request == ki_todo) then
					mpi_ki_selector = .true.
					exit loop_todos
				end if
			end do loop_todos
		end if
		!
		return
	end function


	subroutine print_core_info(n_ki_loc, sum_N_el_loc, min_n_el, max_n_el)
		integer,		intent(in)		::	n_ki_loc
		real(dp),		intent(in)		::	sum_N_el_loc, min_n_el, max_n_el
		!
		write(*,'(a,i3,a,a,a,i8,a)')					"[#",mpi_id,"; core_worker/",cTIME(time()),		&
													"]: ...finished interpolating ",n_ki_loc," kpts"
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		write(*,'(a,i3,a,f4.2,a,f4.2,a,f4.2,a)')	"[#",mpi_id,"; core_worker]: avg el count ",						&
																					sum_N_el_loc/real(n_ki_loc,dp),		&
																				"	(min: ", 							&
																					min_n_el,							&
																				";   max: ", 							&
																					max_n_el,							&
																				")"
		!
		return
	end subroutine







end module core


