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
									use_kspace_ham,			&
									kspace_ham_id,			&
									wf_centers,				&
									use_cart_velo,			&
									kubo_tol, a_latt
	use model_hams,		only:		get_kspace_ham
	use file_io,		only:		write_eig_binary,		&
									write_ham_binary,		&
									write_velo
	use mpi_community,	only:		mpi_root_id, mpi_id, mpi_nProcs, ierr
	use mpi
	use k_space,		only:		get_cart_kpt
	use omp_lib



	implicit none


	private
	public					::		get_wann_interp,		&
									get_kubo_curv	


	contains









!public:
	subroutine get_wann_interp(		do_gauge_trafo, 									&
									H_real, r_real, 									&
									R_frac,	wf_centers,									&
									kpt_idx, kpt, 									&
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
		complex(dp),	allocatable,	intent(in)				::	H_real(:,:,:)
		complex(dp),	allocatable, 	intent(inout)			::	r_real(:,:,:,:)
		integer,						intent(in)				::	kpt_idx
		real(dp),						intent(in)				::	R_frac(:,:), kpt(3)	
		real(dp),		allocatable,	intent(in)				::	wf_centers(:,:)
		real(dp),						intent(out)				::	e_k(:)
		complex(dp),	allocatable,	intent(inout)			::	V_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout)			::	A_ka(:,:,:), Om_kab(:,:,:,:)
		
		complex(dp),	allocatable								::	U_k(:,:), H_ka(:,:,:)
		!
		!
								allocate(	U_k(			size(e_k,1),		size(e_k,1)		)		)		
!		if(	allocated(V_ka)	)	allocate(	H_ka(	3	,	size(e_k,1),		size(e_k,1)	 	)		)
		allocate(	H_ka(	3	,	size(e_k,1),		size(e_k,1)	 	)		)

		!
		!
		!
		if(	.not. use_kspace_ham )	then
			!
			!	INTERPOLATE REAL SPACE BASIS
			call FT_R_to_k(H_real, r_real,  R_frac, wf_centers, kpt_idx, kpt, U_k,  H_ka, A_ka, Om_kab)
		else
			!
			!	USE MODEL HAM 
			call get_kspace_ham(kspace_ham_id,	kpt_idx, kpt, U_k, H_ka)
		end if

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
			if(debug_mode)	call check_H_gauge_herm(kpt_idx, kpt, A_ka, Om_kab, V_ka)
		end if

		!
		return
	end subroutine




	subroutine	get_kubo_curv(	en_k, V_ka, kubo_curv	)
		!
		!	use Kubo formula to calculate curvature of bands
		!
		!		curv_n = -2 sum_np	IM[		V^_   V^_	]	/	(E_np - E_n)**2
		!
		real(dp),						intent(in)			::	en_k(:)
		complex(dp),					intent(in)			::	V_ka(:,:,:)
		complex(dp),	allocatable,	intent(out)			::	kubo_curv(:,:,:)
		real(dp)											::	curv_nn(3,3),	dEE_nnp
		integer												::	n_wf, n, np,  j
		!
		n_wf		=	size(en_k,1)
		allocate(kubo_curv(3,3,n_wf))
		kubo_curv	=	cmplx(0.0_dp,0.0_dp,dp)
		!
		do n = 1, n_wf
			!
			curv_nn	=	0.0_dp
			do np = 1, n_wf
				if(np==n)					cycle
				!
				dEE_nnp				=	(	en_k(np) - en_k(n)	)**2 
				if(	abs(dEE_nnp)< kubo_tol)	cycle
				!
				do j = 1, 3
					kubo_curv(:,j,n)	=	kubo_curv(:,j,n)	+	aimag(	V_ka(:,n,np) * V_ka(j,np,n)	)	/ dEE_nnp
				end do
			end do
		end do
		!
		kubo_curv	=	-2.0_dp * kubo_curv
		!
		return
	end subroutine


!private:
	subroutine FT_R_to_k(H_real, r_real, R_frac, atom_frac, kpt_idx, kpt, H_k,	H_ka, A_ka, Om_kab)			
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
		real(dp),						intent(in)				::	R_frac(:,:), kpt(3)	
		real(dp),		allocatable,	intent(in)				::	atom_frac(:,:)
		integer,						intent(in)				::	kpt_idx
		complex(dp),					intent(out)				::	H_k(:,:)
		complex(dp),	allocatable,	intent(inout)			::	H_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp),		allocatable								::	delta_at(:,:,:)
		real(dp)												::	dR(3), dR_cart(3),kpt_cart(3),ft_angle, two_pi
		complex(dp)												::	ft_phase
		logical													::	use_pos_op, do_en_grad
		integer    												::	sc, a, b, n_sc, n_wf, n, m 
		!
		n_sc	=	size(R_frac,2)
		n_wf	=	size(H_real,1)
		two_pi	=	2.0_dp 	* pi_dp

		!jobs
		do_en_grad		= allocated(H_ka)
		use_pos_op		= allocated(A_ka) .and. allocated(r_real) .and. allocated(Om_kab)
		!
		!init
						H_k		= 	0.0_dp
		if(do_en_grad)	H_ka	= 	0.0_dp
		if(use_pos_op)	A_ka	= 	0.0_dp
		if(use_pos_op)	Om_kab	= 	0.0_dp
		!
		kpt_cart	=	get_cart_kpt(a_latt, kpt)

		!
		! I		INTERPOLATE HAM & VELO(s)
		!
		if( allocated(atom_frac)) then	
			if(	kpt_idx <= mpi_nProcs	)then 
				write(*,'(a,i7.7,a)')	"[#",mpi_id,";FT_R_to_k]: will use TIGHT-B FT CONVENTION"
			end if
			! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!|
			!	CONSIDER ATOMIC POSITIONS IN FOURIER TRANSFORM															!|
			!	THIS IS THE 	"TIGHT BINDING" 	CONVENTION															!|
			!	---------																								!|
			allocate(	delta_at(size(atom_frac,1),size(atom_frac,2),size(atom_frac,2)))
			do n = 1, n_wf
				do m = 1, n_wf
					delta_at(:,n,m)	=	atom_frac(:,m) - atom_frac(:,n)	
				end do
			end do
			!$OMP PARALLEL DEFAULT(none)								&
			!$OMP PRIVATE( 	dR, dR_cart,  ft_angle, ft_phase, m,n)		&
			!$OMP SHARED(  	use_cart_velo,kpt_cart, n_sc, n_wf,kpt_idx, R_frac, delta_at, a_latt, H_real, kpt, H_k, H_ka, do_en_grad, two_pi)
			!$OMP DO REDUCTION(+: H_k, H_ka)
			do sc = 1, n_sc																								!|
				do m = 1 , n_wf																							!|
					do n = 1, n_wf																						!|
						dR(:)		=	R_frac(:,sc) +	delta_at(:,n,m)													!|	
						!	get ft angle																				!|
						if( use_cart_velo	) then																		!|
							dR_cart		=	intern_to_cart(a_latt, dR)													!|			
							ft_angle	=	 			dot_product(	kpt_cart(:),	dR_cart(:)	)! cartesian coords	!|
						else																							!|
							ft_angle	=	two_pi * 	dot_product(		kpt(:) ,		dR(:)	)! internal coords	!|
						end if																							!|
						!																								!|
						!																								!|
						ft_phase	= 	cmplx(	cos(ft_angle), sin(ft_angle),	dp	)	!unit indepedent				!|		
						!Hamilton operator																				!|
						H_k(n,m)	=			H_k(n,m)		+	ft_phase						* H_real(n,m,sc)	!|	
						!																								!|
						!OPTIONAL energy gradients																		!|
						if( do_en_grad)		then																		!|					
							!																							!|	
							!	get dR_cart																				!|										
							if( .not. use_cart_velo ) then																!|																					!|
								dR_cart(:)	=	two_pi  * dR(:)															!|
							end if																						!|
							!																							!|
							H_ka(:,n,m) 	=	H_ka(:,n,m)		+ ft_phase * cmplx(0.0_dp,dR_cart(:),dp)* H_real(n,m,sc)!|
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
			if(	kpt_idx <= mpi_nProcs	)then 
				write(*,'(a,i7.7,a)')	"[",mpi_id,"FT_R_to_k]: will use WANNIER FT CONVENTION"
			end if
			!
			! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!|
			!	NEGLECT ATOMIC POSITIONS																				!|
			!	THIS IS THE 	"WANNIER"	CONVENTION																	!|
			!	---------		
			!$OMP PARALLEL DEFAULT(none)						&
			!$OMP PRIVATE( dR_cart, ft_angle, ft_phase, a,m,n)	&
			!$OMP SHARED(use_cart_velo, n_sc, n_wf,  R_frac, a_latt, H_real, kpt, H_k, H_ka, do_en_grad, two_pi)
			!$OMP DO REDUCTION(+: H_k, H_ka)																						
			do sc = 1, n_sc																								!|
				ft_angle	=	two_pi * dot_product(	kpt(:),R_frac(:,sc)	)		! internal  coord					!|
				ft_phase	= 	cmplx(	cos(ft_angle), sin(ft_angle), dp)												!|	
				!Hamilton operator																						!|
				H_k(:,:)	= 			H_k(:,:)		+	ft_phase								* H_real(:,:,sc)	!|
				!																										!|
				!OPTIONAL energy gradients																				!|
				if( do_en_grad)		then																				!|
					!	get dR_cart																						!|										
					if( use_cart_velo ) then																			!|
						dR_cart(:)	= 	intern_to_cart( a_latt(:,:), R_frac(:,sc)	)									!|
					else																								!|
						dR_cart(:)	=	two_pi * R_frac(:,sc)															!|
					end if																								!|
					!																									!|
					do a = 1, 3																							!|
						H_ka(a,:,:) 	=		H_ka(a,:,:)		+	ft_phase * i_dp * 	dR_cart(a) * H_real(:,:,sc)		!|
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
			if(	kpt_idx <= mpi_nProcs	)then 
				write(*,'(a,i7.7,a)')	"[",mpi_id,"FT_R_to_k]: will interpolate conn & curv additonally"
			end if
			if( allocated(atom_frac)	) then
				! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!|
				!	CONSIDER ATOMIC POSITIONS IN FOURIER TRANSFORM															!|
				!	THIS IS THE 	"TIGHT BINDING" 	CONVENTION															!|
				!	---------
				!$OMP PARALLEL DEFAULT(none) &					
				!$OMP PRIVATE(m,n, dR, dR_cart, ft_angle, ft_phase, a, b) &	
				!$OMP Shared(use_cart_velo, n_sc, n_wf, A_ka, Om_kab, a_latt, R_frac, delta_at, kpt,  r_real, two_pi)	
				!$OMP DO REDUCTION(+: A_ka, Om_kab)																											
				do sc = 1, n_sc																								!|																				
					do m = 1 , n_wf																							!|
						do n = 1, n_wf																						!|
							dR(:)		=	R_frac(:,sc) +	delta_at(:,n,m)				! internal units					!|
							!																								!|
							ft_angle	=	2.0_dp * pi_dp * dot_product(kpt(:),dR(:))										!|
							ft_phase	= 	cmplx(	cos(ft_angle), sin(ft_angle)	,	dp	)								!|
							!																								!|
							!	get dR_cart																					!|										
							if( use_cart_velo ) then																		!|
								dR_cart(:)	= 	intern_to_cart( a_latt(:,:), dR(:)	)										!|
							else																							!|
								dR_cart(:)	=	two_pi * dR(:)																!|
							end if																							!|
							!																								!|									
							!																								!|	
							do a = 1, 3																						!|									
								!connection																					!|					
								A_ka(a,n,m)			=	A_ka(a,n,m)		+	ft_phase				* r_real(a,n,m,sc)		!|
								!curvature																					!|	
								do b = 1, 3																					!|	
									Om_kab(a,b,n,m)	=	Om_kab(a,b,n,m)	&													!|	
													 		+ 	ft_phase * i_dp * dR_cart(a) 	* r_real(b,n,m,sc)			!|	
									Om_kab(a,b,n,m)	=	Om_kab(a,b,n,m)	& 													!|
															- 	ft_phase * i_dp * dR_cart(b) 	* r_real(a,n,m,sc)			!|	
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
				!$OMP PRIVATE(dR_cart, ft_angle, ft_phase, a, b) &	
				!$OMP Shared(use_cart_velo, n_sc, a_latt, R_frac, r_real, A_ka, Om_kab, kpt, two_pi)	
				!$OMP DO REDUCTION(+: A_ka, Om_kab)
				do sc = 1, n_sc																								!|
					ft_angle	=	two_pi * dot_product(kpt(:),R_frac(:,sc))												!|
					ft_phase	= 	cmplx(	cos(ft_angle), sin(ft_angle),	dp	)											!|
					!	get dR_cart																							!|										
					if( use_cart_velo ) then																				!|
						dR_cart(:)	= 	intern_to_cart( a_latt(:,:), R_frac(:,sc)	)										!|
					else																									!|
						dR_cart(:)	=	two_pi * R_frac(:,sc)																!|
					end if																									!|
					!																										!|
					do a = 1, 3																								!|									
						!connection																							!|					
						A_ka(a,:,:)			=	A_ka(a,:,:)		+	ft_phase		* r_real(a,:,:,sc)						!|
						!curvature																							!|	
						do b = 1, 3																							!|	
							Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)	&															!|	
											 		+ 	ft_phase * i_dp * dR_cart(a) * r_real(b,:,:,sc)						!|	
							Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)	& 															!|
													- 	ft_phase * i_dp * dR_cart(b) * r_real(a,:,:,sc)						!|	
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
			call debug_ft_phase(R_frac, delta_at, kpt, a_latt)	
			call write_ham_binary(kpt_idx,	H_k)
			call check_W_gauge_herm(kpt, H_k, H_ka, A_ka, Om_kab)
		end if	
		!
		return
	end subroutine


	pure function intern_to_cart(  a_latt, vec	) result(vec_cart)
		real(dp),	intent(in)		::		a_latt(3,3), vec(3)
		real(dp)					::		vec_cart(3)
		!integer						::		dim
		!
		!vec_cart	=	0.0_dp
		!do dim = 1, 3
		!	vec_cart(:)	=	vec_cart(:)	+	vec(dim)	*	a_latt(dim,:)	
		!end do
		vec_cart(:)		=	matmul(		a_latt, vec(:)	)		
		!
		return
	end function






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
				do n = m+1, n_wf
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
		complex(dp),	allocatable						::	U_dag(:,:), M_in(:,:), work(:,:)
		integer											::	a, b
		!
		allocate(	U_dag(		size(U_k,1), size(U_k,2)	))
		allocate(	M_in(		size(U_k,1), size(U_k,2)	))
		allocate(	work(		size(U_k,1), size(U_k,2)	))
		!
		U_dag	=	conjg(	transpose(	U_k		))
		!
		!
		do a = 1, 3
			!
			!	VELOCITIES
										M_in				=	H_ka(a,:,:)
										H_ka(a,:,:)			=	unit_rot(	M_in,		U_k, U_dag, 	work)

			!	CONNECTION
			if( allocated(A_ka))then	
										M_in				=	A_ka(a,:,:)
										A_ka(a,:,:)			=	unit_rot(	M_in,		U_k, U_dag, 	work)
			end if
			!	CURVATURE
			if( allocated(Om_kab)	)then
				do b = 1,3 
										M_in				=	Om_kab(a,b,:,:)
										Om_kab(a,b,:,:)		=	unit_rot(	M_in,		U_k, U_dag, 	work)
				end do
			end if
			!
		end do
		!
		!
		return
	end subroutine	

	function unit_rot(	M,	U, U_dag, tmp	)	result(M_rot)
		complex(dp),	intent(in)		::	M(:,:), U(:,:), U_dag(:,:) 
		complex(dp),	intent(inout)	::	tmp(:,:)
		complex(dp),	allocatable		::	M_rot(:,:)
		!
		allocate(	M_rot(	size(M,1),	size(M,2)	))
		!
		tmp		=	blas_matmul(	  M		,	U	)
		M_rot	=	blas_matmul(	U_dag	,  tmp	)
		!
		return
	end function



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
		!$OMP PARALLEL DO DEFAULT(none) &
		!$OMP PRIVATE(n,eDiff_mn) &
		!$OMP SHARED(kubo_tol, e_k, H_Ka, D_ka) 
		do m = 1, size(D_ka,3)
			do n = m+1, size(D_ka,2)
				eDiff_mn			=	e_k(m)	- 	e_k(n)
				!
				if(abs(eDiff_mn) > 	kubo_tol	)	then
					D_ka(1:3,n,m)	=	H_ka(1:3,n,m) / 	eDiff_mn
					D_ka(1:3,m,n)	=	H_ka(1:3,m,n) /	( - eDiff_mn	)
				else
					write(*,'(a)',advance="no")	'[;get_gauge_covar_deriv]: '
					write(*,'(a,i6,a,i6)')		'WARNING degenerate bands detetected n=',n,' m=',m
				end if
			end do
		end do
		!$OMP END PARALLEL DO
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
		!$OMP PARALLEL DO DEFAULT(none)	&
		!$OMP COLLAPSE(2)				&
		!$OMP PRIVATE(n,m)				&
		!$OMP SHARED(A_ka, D_ka)
		do m = 1, size(A_ka,3)
			do n = 1, size(A_ka,2)
				A_ka(:,n,m)	=	A_ka(:,n,m)	+ i_dp	*	D_ka(:,n,m)
			end do
		end do
		!$OMP END PARALLEL DO
		!
		return
	end subroutine


	subroutine curv_gaugeTrafo(D_ka, A_ka, Om_kab)
		!	PRB 74, 195118 (2006)	EQ.(27)
		complex(dp),		intent(in)		::	D_ka(:,:,:), A_ka(:,:,:)
		complex(dp),		intent(inout)	::	Om_kab(:,:,:,:)
		integer								::	a, b
		!
		!$OMP PARALLEL DO DEFAULT(none) 	&
		!$OMP COLLAPSE(2)					&
		!$OMP PRIVATE(a,b)					&
		!$OMP SHARED(Om_kab, D_ka, A_ka)
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
		!$OMP END PARALLEL DO
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
	

	subroutine debug_ft_phase(R_frac, delta_at, kpt, a_latt)
		real(dp),						intent(in)				::	R_frac(:,:), delta_at(:,:,:), kpt(3), a_latt(3,3)
		real(dp)												:: 	dR(3), dR_cart(3), kpt_cart(3), &
																	ft_agl, ft_agl_cart, two_pi 
		complex(dp)												::	R_sum(3), R_sum_cart(3),			&
																	ft_phs, ft_phs_cart,&
																	ft_sum, ft_sum_cart, two_pi_i 
		integer													::	m, n, sc, n_sc , n_wf
		!
		n_wf	=	size(delta_at,2)
		n_sc	=	size(R_frac,2)
		!
		kpt_cart(:)	=	get_cart_kpt(a_latt, kpt)
		!
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		if(	mpi_id == mpi_root_id ) then	
			write(*,*)	"[debug_ft_phase]: now check (n,m)=(",n,",",m,")	#sc=",n_sc	,"kpt=",kpt			
			two_pi		=	2.0_dp * pi_dp
			two_pi_i	= 	two_pi * i_dp
			ft_sum 		=	0.0_dp
			ft_sum_cart	=	0.0_dp
			R_sum		=	0.0_dp
			R_sum_cart	=	0.0_dp
			!
			do m = 1 , n_wf	
					do n = 1, n_wf	
			do sc = 1, n_sc										
					
						dR(:)		=	R_frac(:,sc)	+ delta_at(:,n,m)				! 	int						
						dR_cart(:)	=	intern_to_cart( a_latt, dR )	!	cart
						!
						!	INTERNAL UNITS 
						ft_agl		=	two_pi * dot_product(	kpt(:),	dR(:)	)
						ft_phs		=	cmplx( cos(ft_agl), sin(ft_agl), dp)
						ft_sum 		=	ft_sum 		+	ft_phs
						R_sum(:)	=	R_sum(:)	+	two_pi_i	* dR(:)	* ft_phs
						!	-------------------
						!
						!	CARTESIAN COORDS
						ft_agl_cart =	dot_product( kpt_cart(:), dR_cart(:) )
						ft_phs_cart = 	cmplx(	cos(ft_agl_cart), sin(ft_agl_cart), 	dp)
						ft_sum_cart =	ft_sum_cart + 	ft_phs_cart
						R_sum_cart(:)	=	R_sum_cart(:)	+	dR_cart(:)	* ft_phs_cart
						!	-------------------
						!
						end do
						!
			!if( abs(ft_sum - 1.0_dp) > 1e-4_dp) &
			write(*,*)	"[debug_ft_phase]:	WARNING      ft_sum=",ft_sum,'vs ',ft_sum_cart,"=ft_sum_cart"
		
			write(*,*)	"[debug_ft_phase]:	WARNING      R_sum=",R_sum,'.'
						write(*,*)	"[debug_ft_phase]:	WARNING      R_sum_cart=",R_sum_cart,'.'

					end do
				end do
			
		end if
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		!
		return
	end subroutine



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



	subroutine check_W_gauge_herm(kpt, H_k, H_ka, A_ka, Om_kab)
		real(dp),						intent(in)		::		kpt(3)
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
		write(k_string,'(a,f6.2,a,f6.2,a,f6.2,a)')	') at rel. kpt=( ',kpt(1),', ',kpt(2),', ',kpt(3),')'
		
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


	subroutine check_H_gauge_herm(kpt_idx, kpt, A_ka, Om_kab, V_ka)
		integer,						intent(in)			::	kpt_idx
		real(dp),						intent(in)			::	kpt(3)
		complex(dp),	allocatable,	intent(in)			::	A_ka(:,:,:), Om_kab(:,:,:,:)
		complex(dp),					intent(in)			::	V_ka(:,:,:)
		real(dp)											::	max_err
		character(len=31)									::	k_string
		character(len=24)									::	allo_lst
		logical												::	conn, curv, velo, is_herm
		!
		allo_lst	=	" "
		write(k_string,'(a,f6.2,a,f6.2,a,f6.2,a)')	'( ',kpt(1),', ',kpt(2),', ',kpt(3),') '
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