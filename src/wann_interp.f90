!
!
!	ToDo:
!
!	Introduce FFT
!	this decreases the DFT(discrete fourrier transform)
!	from O(n**2) to O(n log(n))
!
!
module wann_interp
	!	module for Wannier interpolation
	!
	!	the interpolation scheme from 
	!			PRB 74, 195118 (2006) 
	!	was used
	use constants,		only:		dp, fp_acc, i_dp, pi_dp	
	use matrix_math,	only:		zheevd_wrapper, 		&
									zheevx_wrapper,			&
									matrix_comm,			& 
									blas_matmul,			&
									is_equal_mat,			&
									is_herm_mat,			&
									is_skew_herm_mat		
	use input_paras,	only:		debug_mode,				&
									kubo_tol
	use file_io,		only:		write_eig_binary,		&
									write_ham_binary
	use mpi_community,	only:		mpi_root_id, mpi_id, mpi_nProcs, ierr
	use omp_lib



	implicit none


	private
	public					::		get_wann_interp	


	contains









!public:
	subroutine get_wann_interp(		do_gauge_trafo, 									&
									H_real, r_real, 									&
									R_frac, atPos, 										&
									kpt_idx, kpt_frac, 									&
									e_k, V_ka, A_ka, Om_kab								&
							)
		!
		!	interpolates the k-space:
		!			-	H_k	: 		hamiltonian
		!			-	H_ka:		k-space derivative of Ham 			
		!			-	A_ka:		Berry conncection (only if r_real given)	
		!			-	Om_kab:		Berry curvature	
		!the real space basis (H_real (real space Hamiltonian) and optionally r_real (real space postition operator)	)
		!
		!	see 	>>>	PRB 74, 195118 (2006)	<<<			for more details on Wannier interpolation
		!
		logical,						intent(in)				::	do_gauge_trafo
		complex(dp),					intent(in)				::	H_real(:,:,:)
		complex(dp),	allocatable, 	intent(inout)			::	r_real(:,:,:,:)
		integer,						intent(in)				::	kpt_idx
		real(dp),						intent(in)				::	R_frac(:,:), kpt_frac(3)	
		real(dp),		allocatable,	intent(in)				::	atPos(:,:)
		real(dp),						intent(out)				::	e_k(:)
		complex(dp),	allocatable,	intent(inout)			::	V_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout)			::	A_ka(:,:,:), Om_kab(:,:,:,:)
		
		complex(dp),	allocatable								::	U_k(:,:), H_ka(:,:,:)
		!
		!
								allocate(	U_k(			size(H_real,1),		size(H_real,2)		)		)		
		if(	allocated(V_ka)	)	allocate(	H_ka(	3	,	size(H_real,1),		size(H_real,2)	 	)		)
		!
		if(			(kpt_idx <= mpi_nProcs)	& 
			.and. 	allocated(atPos)		&
		)then 
					write(*,*)	"[",mpi_id,"get_wann_interp]: will use atomic positions in FT"
		end if
		!
		!
		!ft onto k-space (W)-gauge
		call FT_R_to_k(H_real, r_real,  R_frac, atPos, kpt_idx, kpt_frac, U_k,  H_ka, A_ka, Om_kab)
		!
		!
		!get energies (H)-gauge
		call zheevd_wrapper(U_k, e_k)
		!
		!
		!rotate back to (H)-gauge
		if (allocated( V_ka)) then
			if( do_gauge_trafo	)	call W_to_H_gaugeTRAFO(kpt_idx, e_k, U_k, H_ka, A_ka, Om_kab)
			call get_velo(e_k, H_ka, A_ka, 	V_ka)
			!	DEBUG
			if(debug_mode)	call check_H_gauge_herm(kpt_idx, kpt_frac, A_ka, Om_kab, V_ka)
		end if
		!
		!
		return
	end subroutine







!private:
	subroutine FT_R_to_k(H_real, r_real, R_frac, atom_frac, kpt_idx, kpt_frac, H_k,	H_ka, A_ka, Om_kab)			
		!	interpolates real space Ham and position matrix to k-space,
		!	according to
		!		PRB 74, 195118 (2006)		EQ.(37)-(40)
		!
		!
		!	->	only the H_real, and H_k have to be allocated
		!	->	all other quantities are only calculated if allocated
		!
		complex(dp),					intent(in)				::	H_real(:,:,:)
		complex(dp),	allocatable, 	intent(inout)			::	r_real(:,:,:,:)
		real(dp),						intent(in)				::	R_frac(:,:), kpt_frac(3)	
		real(dp),		allocatable,	intent(in)				::	atom_frac(:,:)
		integer,						intent(in)				::	kpt_idx
		complex(dp),					intent(out)				::	H_k(:,:)
		complex(dp),	allocatable,	intent(inout)			::	H_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp),		allocatable								::	delta_at(:,:,:)
		real(dp)												::	delta_R(3), ft_angle, two_pi
		complex(dp)												::	ft_phase, two_pi_i
		logical													::	use_pos_op, do_en_grad
		integer    												::	sc, a, b, n_sc, n_wf, n, m 
		!
		n_sc	=	size(R_frac,2)
		n_wf	=	size(H_real,1)
		!jobs
		do_en_grad		= allocated(H_ka)
		use_pos_op		= allocated(A_ka) .and. allocated(r_real) .and. allocated(Om_kab)
		!
		!init
						H_k		= 	0.0_dp
		if(do_en_grad)	H_ka	= 	0.0_dp
		if(use_pos_op)	A_ka	= 	0.0_dp
		if(use_pos_op)	Om_kab	= 	0.0_dp
		two_pi					=	2.0_dp * pi_dp
		two_pi_i				=	cmplx(0.0_dp, two_pi, dp)
		!
		!
		! I		INTERPOLATE HAM & VELO(s)
		!
		if( allocated(atom_frac)) then	
			! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!|
			!	CONSIDER ATOMIC POSITIONS IN FOURIER TRANSFORM															!|
			!	THIS IS THE 	"TIGHT BINDING" 	CONVENTION															!|																!|
			!	---------		
			allocate(	delta_at(3,n_wf,n_wf))
			do n = 1, n_wf
				do m = 1, n_wf
					delta_at(:,n,m)	=	atom_frac(:,m) - atom_frac(:,n)	
				end do
			end do
			!$OMP PARALLEL DEFAULT(none)							&
			!$OMP PRIVATE( delta_R,  ft_angle, ft_phase, a,m,n)		&
			!$OMP SHARED(  H_real, H_k, H_ka, R_frac,  delta_at, kpt_frac, n_sc, n_wf, do_en_grad, two_pi, two_pi_i)
			!$OMP DO REDUCTION(+: H_k, H_ka)
			do sc = 1, n_sc																								!|
				do m = 1 , n_wf																							!|
					do n = 1, n_wf																						!|
						delta_R(:)	=	R_frac(:,sc) +	delta_at(:,n,m)						!	internal units			!|
						!																								!|
						ft_angle	=	two_pi * dot_product(kpt_frac(1:3),	delta_R(1:3))	!factor 2pi comes from		!|
						ft_phase	= 	cmplx(	cos(ft_angle), sin(ft_angle)	,	dp	)	! using internanl units		!|
						!																								!|
						!Hamilton operator																				!|
						H_k(n,m)	=			H_k(n,m)		+	ft_phase						* H_real(n,m,sc)	!|	
						!																								!|
						!OPTIONAL energy gradients																		!|
						if( do_en_grad)		then																		!|
							H_ka(:,n,m) 	=	H_ka(:,n,m)		+	ft_phase * two_pi_i	 * delta_R(:) 	* H_real(n,m,sc)!|
						end if																							!|
					end do																								!|
				end do
			end do
			!$OMP END DO
			!$OMP END PARALLEL
			! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!|
			!
			!
		else
			!
			! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!|
			!	NEGLECT ATOMIC POSITIONS																				!|
			!	THIS IS THE 	"WANNIER"	CONVENTION																	!|
			!	---------		
			!$OMP PARALLEL DEFAULT(none)						&
			!$OMP PRIVATE( delta_R,  ft_angle, ft_phase, a,m,n)		&
			!$OMP SHARED(  H_real, H_k, H_ka, R_frac,  kpt_frac, n_sc, n_wf, do_en_grad, two_pi, two_pi_i)
			!$OMP DO REDUCTION(+: H_k, H_ka)																						
			do sc = 1, n_sc
				ft_angle			=	two_pi * dot_product(kpt_frac(1:3), R_frac(:,sc))								!|	
				ft_phase			= 	cmplx(	cos(ft_angle), sin(ft_angle), dp)										!|	
				!Hamilton operator																						!|
				H_k(:,:)			= 			H_k(:,:)		+	ft_phase						* H_real(:,:,sc)	!|
				!																										!|
				!OPTIONAL energy gradients																				!|
				if( do_en_grad)		then																				!|
					do a = 1, 3																							!|
						H_ka(a,:,:) 	=		H_ka(a,:,:)		+	ft_phase * two_pi_i * 	R_frac(a,sc) * H_real(:,:,sc)!|
					end do																								!|
				end if																									!|
			end do
			!$OMP END DO
			!$OMP END PARALLEL
			! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!|
		end if
		!
		!
		! 
		!
		! II	INTERPOLATE CONN & CURV
		!
		if(use_pos_op) then
			if( allocated(atom_frac)	) then
				! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!|
				!	CONSIDER ATOMIC POSITIONS IN FOURIER TRANSFORM															!|
				!	THIS IS THE 	"TIGHT BINDING" 	CONVENTION															!|
				!	---------
				!$OMP PARALLEL DEFAULT(none) &					
				!$OMP PRIVATE(m,n, delta_R, ft_angle, ft_phase, a, b) &	
				!$OMP Shared(n_sc, n_wf, A_ka, Om_kab, R_frac, delta_at, kpt_frac, r_real, two_pi, two_pi_i)	
				!$OMP DO REDUCTION(+: A_ka, Om_kab)																											
				do sc = 1, n_sc																								!|																				
					do m = 1 , n_wf																							!|
						do n = 1, n_wf																						!|
							delta_R(:)	=	R_frac(:,sc) +	delta_at(:,n,m)								!	internal units	!|
							!																								!|
							ft_angle	=	two_pi * dot_product(kpt_frac(1:3),	delta_R(1:3))								!|
							ft_phase	= 	cmplx(	cos(ft_angle), sin(ft_angle)	,	dp	)								!|
							!																								!|									
							!																								!|	
							do a = 1, 3																						!|									
								!connection																					!|					
								A_ka(a,:,:)			=	A_ka(a,:,:)		+	ft_phase				* r_real(a,:,:,sc)		!|
								!curvature																					!|	
								do b = 1, 3																					!|	
									Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)	&													!|	
													 		+ 	ft_phase * two_pi_i * delta_R(a) 	* r_real(b,:,:,sc)		!|	
									Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)	& 													!|
															- 	ft_phase * two_pi_i * delta_R(b) 	* r_real(a,:,:,sc)		!|	
								end do																						!|
							end do																							!|
						end do																								!|
					end do																									!|
					!																										!|
					!																							!~~~~~~~~~~~!|
				end do 					
				!$OMP END DO
				!$OMP END PARALLEL 
				! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!|
				!
				!
			else
				! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!|
				!	NEGLECT ATOMIC POSITIONS																				!|
				!	THIS IS THE 	"WANNIER"	CONVENTION																	!|
				!	---------
				!$OMP PARALLEL DEFAULT(none) &					
				!$OMP PRIVATE(ft_angle, ft_phase, a, b) &	
				!$OMP Shared(n_sc, A_ka, Om_kab, R_frac, kpt_frac, r_real, two_pi, two_pi_i)	
				!$OMP DO REDUCTION(+: A_ka, Om_kab)
				do sc = 1, n_sc																								!|
					ft_angle	=	two_pi * dot_product(kpt_frac(1:3),	R_frac(1:3,sc))										!|
					ft_phase	= 	cmplx(	cos(ft_angle), sin(ft_angle)	,	dp	)										!|
					!																										!|									
					!																										!|	
					do a = 1, 3																								!|									
						!connection																							!|					
						A_ka(a,:,:)			=	A_ka(a,:,:)		+	ft_phase					* r_real(a,:,:,sc)			!|
						!curvature																							!|	
						do b = 1, 3																							!|	
							Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)	&															!|	
											 		+ 	ft_phase * two_pi_i * R_frac(a,sc) * r_real(b,:,:,sc)				!|	
							Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)	& 															!|
													- 	ft_phase * two_pi_i * R_frac(b,sc) * r_real(a,:,:,sc)				!|	
						end do																								!|
					end do																									!|
					!																										!|
					!																							!~~~~~~~~~~~!|
				end do 																										
				!$OMP END DO
				!$OMP END PARALLEL
				! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!|
			end if
		end if	
		!
		if(debug_mode) then	
			call write_ham_binary(kpt_idx,	H_k)
			call check_W_gauge_herm(kpt_frac, H_k, H_ka, A_ka, Om_kab)
		end if	
		!
		return
	end subroutine







	subroutine get_velo(en_k, H_ka, A_ka, V_ka)
		!	calc the (H)-gauge velocity matrix
		!
		!	PRB 74, 195118 (2006) EQ.(31)
		real(dp),						intent(in)		::		en_k(:)
		complex(dp),					intent(in)		::		H_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout) 	::		A_ka(:,:,:)
		complex(dp),					intent(out)		::		V_ka(:,:,:)
		complex(dp)										::		cmplx_eDiff
		integer											::		m, n, n_wf
		!
		n_wf		=	size(en_k,1)
		!
		!
		if( allocated(A_ka)	) then
			!	if CONNECTION given 
			!	-> apply correction to H_Ka
			!
			!$OMP PARALLEL DO DEFAULT(none) &
			!$OMP PRIVATE(n,cmplx_eDiff) &
			!$OMP SHARED(n_wf, en_k, H_Ka, A_ka, V_Ka) 
			do m = 1, n_wf
				do n = 1, n_wf
					if(	n >	m )	then
						cmplx_eDiff	=	cmplx(	0.0_dp	,	en_k(m) - en_k(n),	dp)
						!
						V_ka(:,n,m)	= H_ka(:,n,m)	- cmplx_eDiff	* 	A_ka(:,n,m)
						V_ka(:,m,n)	= H_ka(:,m,n)	+ cmplx_eDiff	*	A_ka(:,m,n)
					end if
				end do
			end do
			!$OMP END PARALLEL DO
		else
			V_ka		=	H_ka
		end if
		!
		!
		return
	end subroutine



!
!
!
!
!
!
!
!	**************************************************************************************************************************************************
!	**************************************************************************************************************************************************
!	**************************************************************************************************************************************************
!
!						~~~~			GAUGE TRAFO ROUTINES				~~~~~~~~~~~~~~
!
!	**************************************************************************************************************************************************
!	**************************************************************************************************************************************************
!	**************************************************************************************************************************************************
!
	subroutine W_to_H_gaugeTRAFO(kpt_idx, e_k, U_k, H_ka, A_ka, Om_kab)
		!	see PRB 74, 195118 (2006) EQ. 21 - 31
		integer,						intent(in)		::	kpt_idx
		real(dp),						intent(in)		::	e_k(:) 
		complex(dp),					intent(in)		::	U_k(:,:)
		complex(dp),	allocatable,	intent(inout)	::	H_ka(:,:,:)
		complex(dp),	allocatable, 	intent(inout)	::	A_ka(:,:,:), Om_kab(:,:,:,:)
		complex(dp),	allocatable						::	D_ka(:,:,:) 
		!
		!
		!debug (W)-gauge
		if(debug_mode)	then
			call check_velo(U_k, H_ka)
			call write_eig_binary(kpt_idx,	U_k)
		end if
		!
		!
		!do 	(W) -> (Hbar)
		call rotate_gauge(U_k, H_ka, 	A_ka, Om_kab )
		if(debug_mode)	call check_Hbar_gauge_herm(H_ka, A_ka, Om_kab)
		!
		!
		!conn/curv	 (Hbar) -> (H)
		if( allocated(A_ka) )	then
			!	need gauge covar deriv (remember A_ka is not gauge covariant)
			allocate(		D_ka(	3,	size(H_ka,2),	size(H_ka,3))				)
			!
			!
			call get_gauge_covar_deriv(e_k, H_ka, D_ka)
			!
			call conn_gaugeTrafo(D_ka, A_ka)
			call curv_gaugeTrafo(D_ka, A_ka, Om_kab)
		end if
		!
		!
	end subroutine

	subroutine rotate_gauge(U_k, H_ka, A_ka, Om_kab)
		!	PRB 74, 195118 (2006)	EQ.(21)
		complex(dp),					intent(in)		::	U_k(:,:)
		complex(dp), 					intent(inout)	::	H_ka(:,:,:)
		complex(dp), 	allocatable,	intent(inout)	::	A_ka(:,:,:), Om_kab(:,:,:,:)
		complex(dp),	allocatable						::	U_dag(:,:), tmp(:,:), M_in(:,:)
		integer											::	a, b
		!
		allocate(	U_dag(		size(U_k,1), size(U_k,2)	))
		allocate(	tmp(		size(U_k,1), size(U_k,2)	))
		allocate(	M_in(		size(U_k,1), size(U_k,2)	))

		U_dag	=	conjg(	transpose(	U_k		))
		!
		!
		do a = 1, 3
			!
			!	VELOCITIES
										M_in			=	H_ka(a,:,:)
										call	blas_matmul(	U_dag, M_in, 		tmp			)
										call	blas_matmul(	tmp, U_k,			M_in		)
										H_ka(a,:,:)		=	M_in(:,:)
			!	CONNECTION
			if( allocated(A_ka))then	
										tmp				=	A_ka(a,:,:)
										tmp(:,:)		=	matmul(	U_dag		,	tmp		)
										A_ka(a,:,:)		=	matmul(	tmp			, 	U_k		)
			end if
			!	CURVATURE
			if( allocated(Om_kab)	)then
				do b = 1,3 
										tmp				=	Om_kab(a,b,:,:)
										tmp				=	matmul(	U_dag		,	tmp		)
										Om_kab(a,b,:,:)	=	matmul(	tmp			,	U_k		)	
				end do
			end if
			!
		end do
		!
		!
		return
	end subroutine	



	subroutine get_gauge_covar_deriv(e_k, H_ka,	D_ka )
		!	PRB 74, 195118 (2006)	EQ.(24)
		real(dp),			intent(in)		::	e_k(:)
		complex(dp),		intent(in)		::	H_ka(:,:,:)
		complex(dp),		intent(out)		::	D_ka(:,:,:)
		integer								::	m, n, a
		real(dp)							::	eDiff_mn, max_err
		!
		D_ka(:,:,:)	=	cmplx(0.0_dp, 0.0_dp, dp)
		!
		do m = 1, size(D_ka,3)
			do n = 1, size(D_ka,2)
				if(	n >	m )	then
					!
					!
					eDiff_mn	=	e_k(m)	- e_k(n)
					!
					if(abs(eDiff_mn) > 	kubo_tol	)	then
						D_ka(1:3,n,m)	=	H_ka(1:3,n,m) / 	eDiff_mn
						D_ka(1:3,m,n)	=	H_ka(1:3,m,n) /	( - eDiff_mn	)
					else
						write(*,'(a)',advance="no")	'[;get_gauge_covar_deriv]: '
						write(*,'(a,i6,a,i6)')		'WARNING degenerate bands detetected n=',n,' m=',m
					end if
					!
					!
				end if
			end do
		end do
		!
		if(debug_mode)	then
			do a = 1, 3
				!
				!
				if( .not. is_skew_herm_mat(D_ka(a,:,:), max_err)	)	then
					write(*,'(a,i1,a,f16.7)')	"[get_gauge_covar_deriv/DEBUG-MODE]:	WARNING D_(k,a=",a,&
																						") is not skew hermitian, max_err=", max_err
					if(.not. is_herm_mat(H_ka(a,:,:),max_err)	) then
						write(*,'(a,i1,a,f16.7)')	"[get_gauge_covar_deriv/DEBUG-MODE]:	WARNING H_(k,a=",a,&
																						") is not hermitian, max_err=", max_err
					end if
				end if
				!
				!
			end do
		end if
		!
		return
	end subroutine


	pure subroutine conn_gaugeTrafo(D_ka, A_ka)
		!	PRB 74, 195118 (2006)	EQ.(25)
		!
		!	Lapack
		!		https://software.intel.com/en-us/mkl-developer-reference-fortran-gemm#90EAA001-D4C8-4211-9EA0-B62F5ADE9CF0
		!		C :- 
		complex(dp),		intent(in)		::	D_ka(:,:,:)
		complex(dp),		intent(inout)	::	A_ka(:,:,:)
		!
		!
		A_ka	=	A_ka	+ i_dp	*	D_ka
		!
		return
	end subroutine


	subroutine curv_gaugeTrafo(D_ka, A_ka, Om_kab)
		!	PRB 74, 195118 (2006)	EQ.(27)
		complex(dp),		intent(in)		::	D_ka(:,:,:), A_ka(:,:,:)
		complex(dp),		intent(inout)	::	Om_kab(:,:,:,:)
		complex(dp),	allocatable			::	mat_comm(:,:)
		integer								::	a, b
		!
		allocate(	mat_comm(	size(Om_kab,3),size(Om_kab,4)	)		)
		!
		do b = 1, 3
			do a = 1, 3
				Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)		-			matrix_comm(	D_ka(a,:,:), 	A_ka(b,:,:)		)
				!
				!
				Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)		+			matrix_comm(	D_ka(b,:,:), 	A_ka(a,:,:)		)
				!
				Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)		-	i_dp *	matrix_comm( D_ka(a,:,:), 	D_ka(b,:,:))
			end do
		end do
		!
		!
		return 
	end subroutine











!
!
!
!
!
!
!
!	**************************************************************************************************************************************************
!	**************************************************************************************************************************************************
!	**************************************************************************************************************************************************
!
!						~~~~			DEBUGGING HELPER ROUTINES				~~~~~~~~~~~~~~
!
!	**************************************************************************************************************************************************
!	**************************************************************************************************************************************************
!	**************************************************************************************************************************************************
!



	logical function curv_is_herm( Om_kab, max_err)
		complex(dp), allocatable,	intent(in)		::	Om_kab(:,:,:,:)
		real(dp),			intent(out)		::	max_err
		integer								::	a, b
		!
		!
		if(allocated(Om_kab))then 
			if( size(Om_kab,3) == size(Om_kab,4) )	then
				curv_is_herm	=	.true.
				do a = 1, 3
					do b = 1, 3
						if(.not. is_herm_mat(Om_kab(b,a,:,:),max_err)) 	then
							curv_is_herm = .false.
							!write(*,*)	"[conn_curv_is_herm/DEBUG-MODE]: Om_k,a=",b,",b=",a," is not hermitian"
						end if
					end do
				end do
			else
				curv_is_herm	= .false.
				stop " [curv_is_herm/DEBUG-MODE]: k-space matrices life on different basis sets"
			end if
		end if
		!
		return
	end function

	logical function velo_is_herm(V_ka, max_err)
		complex(dp),		intent(in)		::	V_ka(:,:,:)
		real(dp),			intent(out)		::	max_err
		integer								::	a
		!
		velo_is_herm	=	 ( size(V_ka,2) == size(V_ka,3) 	) .and. (size(V_ka,1)==3)
		!
		!
		if(velo_is_herm)	then
			do a = 1, 3
				if(.not. is_herm_mat(V_ka(a,:,:),max_err)) 	then
					velo_is_herm = .false.
				end if
			end do
		else
			velo_is_herm	= .false.
			stop " [velo_is_herm/DEBUG-MODE]: ERROR velo operator matrix not symmetric"
		end if
		!
		return
	end function


	subroutine check_velo(U_k, VW_ka)
		!	check if rotating forward and backwards wants is identity operation
		!
		complex(dp),		intent(in)		::	U_k(:,:), VW_ka(:,:,:)
		complex(dp),	allocatable			::	VW_new_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:)
		!
		allocate(		VW_new_ka(size(VW_ka,1), size(VW_ka,2),size(VW_ka,3))	)
		!
		VW_new_ka	=	VW_ka
		!
		call rotate_gauge(U_k, VW_new_ka, A_ka, Om_kab)
		!
		!	now rotate back to wannier gauge
		call rotate_gauge(	conjg(transpose(U_k)),	VW_new_ka,	A_ka, Om_kab)
		!
		!
		if(			is_equal_mat(	1e-9_dp	, VW_ka(1,:,:), VW_new_ka(1,:,:))	&
			.and.	is_equal_mat(	1e-9_dp	, VW_ka(2,:,:), VW_new_ka(2,:,:))	&
			.and.	is_equal_mat(	1e-9_dp	, VW_ka(3,:,:), VW_new_ka(3,:,:))	&
			)	then
			write(*,*)	"[wann_interp/check_velo]:	SUCCESS gauge consistency seems fine "
		else
			write(*,*)	"[wann_interp/check_velo]:	WARNING gauge consistency not given "
		end if

		return
	end subroutine



	subroutine check_W_gauge_herm(kpt_frac, H_k, H_ka, A_ka, Om_kab)
		real(dp),						intent(in)		::		kpt_frac(3)
		complex(dp),					intent(in)		::		H_k(:,:)
		complex(dp),	allocatable, 	intent(in)		::		H_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:) 
		real(dp)										::		max_err	
		character(len=44)								::		k_string
		character(len=51)								::		warn_msg
		character(len=26)								::		max_string
		character(len=32)								::		allo_lst
		logical 										::		ham_herm, velo_herm, conn_herm, curv_herm
		!
		warn_msg	=								'[check_w_gauge_herm/DEBUG-MODE]:	WARNING (W)-gauge '
		max_string	=								' IS NOT hermitian(max_err='
		write(k_string,'(a,f6.2,a,f6.2,a,f6.2,a)')	') at rel. kpt=( ',kpt_frac(1),', ',kpt_frac(2),', ',kpt_frac(3),')'
		
		!
		!	CHECK HAMILTONIAN
		ham_herm	=	 is_herm_mat( H_k, max_err)
		if(.not. ham_herm) 				write(*,'(a,e16.7,a)')	warn_msg//"H_k"//max_string,	max_err ,k_string 
		allo_lst	=	" ham, "
		!
		!	CHECK DERIVATIVE OF HAM
		if(allocated(H_ka)) then
			allo_lst	=	trim(allo_lst)	//	"velo, "
			velo_herm 	=	 velo_is_herm( H_ka, max_err) 
			if(.not. velo_herm) 		write(*,'(a,e16.7,a)')	warn_msg//"H_ka"//max_string,	max_err ,k_string
		else
			velo_herm	=	.true.
		end if
		!
		!	CHECK CONNECTION
		if(allocated(A_ka)) then
			allo_lst	=	trim(allo_lst)	//	"conn, "
			conn_herm	=	velo_is_herm( A_ka,	max_err)
			if(.not. conn_herm)			write(*,'(a,e16.7,a)')	warn_msg//"A_ka"//max_string,	max_err ,k_string
		else
			conn_herm = .true.
		end if
		!
		!	CHECK CURVATURE
		if(allocated(Om_kab))	then
			allo_lst	=	trim(allo_lst)	//	"curv"
			curv_herm	= curv_is_herm( Om_kab,	max_err)
			if (.not. curv_herm) 		write(*,'(a,e16.7,a)')	warn_msg//"Om_kab"//max_string,	max_err ,k_string
		else
			curv_herm = .true.
		end if
		!
		!
		if( (ham_herm .and. velo_herm) .and. (conn_herm .and. curv_herm) ) then
			write(*,*)	"[check_w_gauge_herm/DEBUG-MODE]: SUCCESS (W)-gauge quantities ("//trim(allo_lst)//") are hermitian"
		end if
		!
		return
	end subroutine


	subroutine check_Hbar_gauge_herm(H_ka, A_ka, Om_kab)
		complex(dp),				intent(in)			::	H_ka(:,:,:) 
		complex(dp), allocatable,	intent(in)			::	A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp)										::	max_err
		!
		if( .not. velo_is_herm(H_ka,max_err))	then
			write(*,'(a,f16.7)')	'[check_Hbar_gauge_herm/DEBUG-MODE]: WARNING	(Hbar)-gauge H_ka IS NOT hermitian, max_err=',max_err	
		end if
		!
		if(allocated(A_ka)) then
			if( .not. velo_is_herm(A_ka,max_err))	then
				write(*,'(a,f16.7)')	'[check_Hbar_gauge_herm/DEBUG-MODE]: WARNING	(Hbar)-gauge A_ka IS NOT hermitian, max_err=',max_err	
			end if
		end if
		!
		if(allocated(Om_kab))	then
			if( .not. curv_is_herm(Om_kab,max_err) )	then
				write(*,'(a,f16.7)')	'[check_Hbar_gauge_herm/DEBUG-MODE]: WARNING	(Hbar)-gauge Om_kab IS NOT hermitian, max_err=',max_err				
			end if
		end if
		!
		!
		return
	end subroutine


	subroutine check_H_gauge_herm(kpt_idx, kpt_frac, A_ka, Om_kab, V_ka)
		integer,						intent(in)			::	kpt_idx
		real(dp),						intent(in)			::	kpt_frac(3)
		complex(dp),	allocatable,	intent(in)			::	A_ka(:,:,:), Om_kab(:,:,:,:)
		complex(dp),					intent(in)			::	V_ka(:,:,:)
		real(dp)											::	max_err
		character(len=31)									::	k_string
		character(len=24)									::	allo_lst
		logical												::	conn, curv, velo, is_herm
		!
		allo_lst	=	" "
		write(k_string,'(a,f6.2,a,f6.2,a,f6.2,a)')	'( ',kpt_frac(1),', ',kpt_frac(2),', ',kpt_frac(3),') '
		is_herm	= .true.
		!
		!	VELOCITY
		allo_lst	=	trim(allo_lst) // "velo"
		velo		=	velo_is_herm( V_ka, max_err )
		is_herm 	=	is_herm .and. velo
		if(.not. velo) 		write(*,'(a,f16.7)')	"[check_H_gauge_herm/DEBUG-MODE]:	"								//	&
														"WARNING (H)-gauge V_KA IS NOT hermitian at rel. kpt= "			//	&
														k_string//"max_err=", max_err
		!
		!	CONNECTION
		if(allocated(A_ka)) then
			allo_lst	=	trim(allo_lst) // "conn, "
			conn	=	velo_is_herm(A_ka, max_err)
			is_herm =	conn
			if(.not. conn) 		write(*,'(a,f16.7)')	"[check_H_gauge_herm/DEBUG-MODE]:	"								//	&
															"WARNING (H)-gauge A_ka IS NOT hermitian at rel. kpt= "			//	&
															k_string//"max_err=", max_err
		else
			write(*,'(a,f16.7)')	"[check_H_gauge_herm/DEBUG-MODE]: NOTE	connection was not calculated"
		end if
		!
		!
		!	CURVATURE
		if(allocated(Om_kab)) then
			allo_lst	=	trim(allo_lst) // "curv, "
			curv	=	curv_is_herm( Om_kab, max_err)
			is_herm =	is_herm .and. curv
			if(.not. curv) 		write(*,'(a,f16.7)')	"[check_H_gauge_herm/DEBUG-MODE]:	"						 		//	&
															"WARNING (H)-gauge Om_kab IS NOT hermitian at rel. kpt= "		//	&
															k_string//"max_err=", max_err
		else
			write(*,'(a,f16.7)')	"[check_H_gauge_herm/DEBUG-MODE]: NOTE	curvature was not calculated"
		end if
		!
		!
		!	SUCCESS MESSAGE
		if(	is_herm ) 	write(*,'(a,i8)')	"[check_H_gauge_herm/DEBUG-MODE]:	"						//	&
															"SUCCESS (H)-gauge quantities ("//trim(allo_lst)//") are hermitian "	//  &
															" at  #kpt= ", kpt_idx		
		!
		return
	end subroutine










end module wann_interp