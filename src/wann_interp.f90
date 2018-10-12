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
	use constants,		only:		dp, fp_acc, i_dp	
	use matrix_math,	only:		zheevd_wrapper, 		&
									matrix_comm,			& 
									uni_gauge_trafo,		&

									is_herm_mat,			&
									is_skew_herm_mat		
	use input_paras,	only:		kubo_tol, debug_mode






	implicit none


	private
	public					::		get_wann_interp	


	contains



!public:
	subroutine get_wann_interp(do_gauge_trafo, H_real, r_real, a_latt, recip_latt, R_frac, kpt_rel, 	e_k, V_ka, A_ka, Om_kab  )
		logical,						intent(in)				::	do_gauge_trafo
		complex(dp),					intent(in)				::	H_real(:,:,:)
		complex(dp),	allocatable, 	intent(inout)			::	r_real(:,:,:,:)
		real(dp),						intent(in)				::	a_latt(3,3), recip_latt(3,3),	& 
																	R_frac(:,:), kpt_rel(3)	
		real(dp),						intent(out)				::	e_k(:)
		complex(dp),	allocatable,	intent(inout)			::	V_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout)			::	A_ka(:,:,:), Om_kab(:,:,:,:)
		
		complex(dp),	allocatable								::	U_k(:,:), H_ka(:,:,:)
		!
		!
								allocate(	U_k(			size(H_real,1),		size(H_real,2)		)		)		
		if(	allocated(V_ka)	)	allocate(	H_ka(	3	,	size(H_real,1),		size(H_real,2)	 	)		)
		!
		!
		!ft onto k-space (W)-gauge
		call FT_R_to_k(H_real, r_real, a_latt, recip_latt, R_frac, kpt_rel, U_k,  H_ka, A_ka, Om_kab)
		!
		!get energies (H)-gauge
		call zheevd_wrapper(U_k, e_k)
		!
		!rotate back to (H)-gauge
		if( allocated(V_ka)	.and. do_gauge_trafo	)			call W_to_H_gaugeTRAFO(e_k, U_k, H_ka, A_ka, Om_kab)
		if( 				.not. do_gauge_trafo	)			write(*,*)	"[get_wann_interp]: WARNING, Wannier gauge (W) selected!"
		if(	allocated(V_ka)							)			call get_velo(e_k, H_ka, A_ka, 	V_ka)
		!
		!
		!	DEBUG
		if(debug_mode)	call check_H_gauge_herm(kpt_rel, H_ka, A_ka, Om_kab, V_ka)
		!
		return
	end subroutine







!private
	subroutine FT_R_to_k(H_real, r_real, a_latt, recip_latt, R_frac, kpt_rel, H_k,	H_ka, A_ka, Om_kab)			
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
		real(dp),						intent(in)				::	a_latt(3,3), recip_latt(3,3),	&
																	R_frac(:,:), kpt_rel(3)	
		complex(dp),					intent(out)				::	H_k(:,:)
		complex(dp),	allocatable,	intent(inout)			::	H_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp)												::	R_abs(3), kpt_abs(3), ft_angle
		complex(dp)												::	ft_phase
		logical													::	use_pos_op, do_en_grad
		integer    												::	sc, a, b 
		!
		!jobs
		do_en_grad		= allocated(H_ka)
		use_pos_op		= allocated(A_ka) .and. allocated(r_real) .and. allocated(Om_kab)
		!
		!init
		kpt_abs		= 	matmul(	recip_latt	, kpt_rel	)	
						H_k		= 0.0_dp
		if(do_en_grad)	H_ka	= 0.0_dp
		if(use_pos_op)	A_ka	= 0.0_dp
		if(use_pos_op)	Om_kab	= 0.0_dp			
		!
		!sum real space cells
		do sc = 1, size(R_frac,2)
			R_abs(:)	=	matmul(	a_latt(:,:),	R_frac(:,sc) )
			ft_angle	=	dot_product(kpt_abs(1:3),	R_abs(1:3))
			ft_phase	= 	cmplx(	cos(ft_angle), sin(ft_angle)	,	dp	)
			!
			ft_phase	=	cmplx(1.0_dp, 0.0_dp, dp)
			!
			!Hamilton operator
			H_k(:,:)				= 	H_k(:,:)			+	ft_phase 					* H_real(:,:,sc)	
			!
			do a = 1, 3
				!OPTIONAL energy gradients
				if( do_en_grad)		then
					H_ka(a,:,:) 		=	H_ka(a,:,:)		+	ft_phase * i_dp * R_abs(a) * H_real(:,:,sc)
				end if
				!OPTIONAL position operator
				if( use_pos_op )	then
					!connection
					A_ka(a,:,:)			=	A_ka(a,:,:)		+	ft_phase					* r_real(a,:,:,sc)
					!curvature
					do b = 1, 3
						Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:) + 	ft_phase * i_dp * R_abs(a) * r_real(b,:,:,sc)
						Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:) - 	ft_phase * i_dp * R_abs(b) * r_real(a,:,:,sc)
					end do
				end if 
			end do
			!
			!
		end do		
		!
		!
		if(debug_mode) 	call check_W_gauge_herm(kpt_rel, H_k, H_Ka, A_ka, Om_kab)
		!
		return
	end subroutine





!gauge TRAFOs
	subroutine W_to_H_gaugeTRAFO(e_k, U_k, H_ka, A_ka, Om_kab)
		!	see PRB 74, 195118 (2006) EQ. 21 - 31
		real(dp),						intent(in)		::	e_k(:) 
		complex(dp),					intent(in)		::	U_k(:,:)
		complex(dp),					intent(inout)	::	H_ka(:,:,:)
		complex(dp),	allocatable, 	intent(inout)	::	A_ka(:,:,:), Om_kab(:,:,:,:)
		complex(dp),	allocatable						::	D_ka(:,:,:) 
		!
		!do 	(W) -> (Hbar)
		call rotate_gauge(U_k, H_ka, 	A_ka, Om_kab )
		!
		if(debug_mode)	call check_Hbar_gauge_herm(A_ka, Om_kab)
		!
		!conn/curv	 (Hbar) -> (H)
		if( allocated(A_ka) )	then
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
		integer											::	a, b
		!
		do a = 1, 3
										call uni_gauge_trafo( U_k,		H_ka(a,:,:)		)
			if( allocated(A_ka)		)	call uni_gauge_trafo( U_k, 		A_ka(a,:,:)		)
			if( allocated(Om_kab)	)then
				do b = 1,3 
										call uni_gauge_trafo( U_k,		Om_kab(a,b,:,:)	)
				end do
			end if
		end do
		!
		return
	end subroutine	


	pure subroutine get_velo(e_k, H_ka, A_ka, V_k)
		!	calc the (H)-gauge velocity matrix
		!
		!	PRB 74, 195118 (2006) EQ.(31)
		real(dp),						intent(in)		::		e_k(:)
		complex(dp),					intent(in)		::		H_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout) 	::		A_ka(:,:,:)
		complex(dp),					intent(out)		::		V_k(:,:,:)
		complex(dp)										::		eDiff
		integer											::		m, n
		!
		V_k		=	H_ka
		!
		!
		if( allocated(A_ka)	) then
			do m = 1, size(V_k,3)
				do n = 1, size(V_k,2)
					if(	n/=	m)	then
						eDiff	=	cmplx(		e_k(m) - e_k(n),		0.0_dp,	dp)
						!
						V_k(:,n,m)	= V_k(:,n,m)	- i_dp	* 	eDiff	* 	A_ka(:,n,m)
					end if
				end do
			end do
		end if
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
		real(dp)							::	eDiff_mn
		!
		D_ka(:,:,:)	=	cmplx(0.0_dp, 0.0_dp, dp)
		!
		do m = 1, size(D_ka,3)
			do n = 1, size(D_ka,2)
				if(	n/=	m )	then
					!
					!
					eDiff_mn	=	e_k(m)	- e_k(n)
					!
					if(abs(eDiff_mn) < 	kubo_tol	)	then
						eDiff_mn	= sign(kubo_tol,eDiff_mn)
						write(*,'(a)',advance="no")	'[;get_gauge_covar_deriv]: '
						write(*,'(a,i6,a,i6)')		'WARNING degenerate bands detetected n=',n,' m=',m
					end if
					!
					!
					D_ka(1:3,n,m)	=	H_ka(1:3,n,m) / 	eDiff_mn
				end if
			end do
		end do
		!
		if(debug_mode)	then
			do a = 1, 3
				if( .not. is_skew_herm_mat(D_ka(a,:,:))		 )	then
					write(*,'(a,i1,a)')	"[get_gauge_covar_deriv/DEBUG-MODE]:	WARNING D_k,a=",a," is not skew hermitian"
				end if
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



!DEBUG
	logical function curv_is_herm( Om_kab)
		complex(dp),		intent(in)		::	Om_kab(:,:,:,:)
		integer								::	a, b
		!
		!
		if(size(Om_kab,3) == size(Om_kab,4) )	then
			curv_is_herm	=	.true.
			do a = 1, 3
				do b = 1, 3
					if(.not. is_herm_mat(Om_kab(b,a,:,:))) 	then
						curv_is_herm = .false.
						!write(*,*)	"[conn_curv_is_herm/DEBUG-MODE]: Om_k,a=",b,",b=",a," is not hermitian"
					end if
				end do
			end do
		else
			curv_is_herm	= .false.
			stop " [curv_is_herm/DEBUG-MODE]: k-space matrices life on different basis sets"
		end if
		!
		return
	end function

		logical function velo_is_herm(V_ka)
		complex(dp),		intent(in)		::	V_ka(:,:,:)
		integer								::	a
		logical								::	basis, velo
		!
		basis	=	size(V_ka,2) == size(V_ka,3) 	
		!
		!
		if(basis)	then
			velo_is_herm	=	.true.
			do a = 1, 3
				if(.not. is_herm_mat(V_ka(a,:,:))) 	then
					velo = .false.
				end if
			end do
		else
			velo_is_herm	= .false.
			stop " [velo_is_herm/DEBUG-MODE]: ERROR velo operator matrix not symmetric"
		end if
		!
		return
	end function



	subroutine check_W_gauge_herm(kpt_rel, H_k, H_ka, A_ka, Om_kab)
		real(dp),			intent(in)			::		kpt_rel(3)	
		complex(dp),		intent(in)			::		H_k(:,:), H_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:) 
		character(len=30)					::	k_string
		!
		write(k_string,'(a,f6.2,a,f6.2,a,f6.2,a)')	'( ',kpt_rel(1),', ',kpt_rel(2),', ',kpt_rel(3),')'
		!
		if(	.not. is_herm_mat(H_k) ) then
			write(*,'(a,200(f6.2))')	"[check_w_gauge_herm/DEBUG-MODE]:	WARNING (W)-gauge H_k IS NOT hermitian at rel. kpt="//k_string
		end if
		!
		if( .not. velo_is_herm(H_ka)) then
			write(*,'(a,200(f6.2))')	"[check_w_gauge_herm/DEBUG-MODE]:	WARNING (W)-gauge H_ka IS NOT hermitian at rel. kpt="//k_string
		end if
		!
		if( .not. velo_is_herm( A_ka)	) 	then
			write(*,*)					"[check_w_gauge_herm/DEBUG-MODE]:	WARNING (W)-gauge A_ka IS NOT hermitian at rel. kpt= "//k_string
		end if
		!
		if( .not. curv_is_herm( Om_kab)	)	then
			write(*,*)					"[check_w_gauge_herm/DEBUG-MODE]:	WARNING (W)-gauge OM_kab IS NOT hermitian at rel. kpt= "//k_string
		end if
		!
		!
		return
	end subroutine

	subroutine check_Hbar_gauge_herm(A_ka, Om_kab)
		complex(dp),	intent(in)			::	 A_ka(:,:,:), Om_kab(:,:,:,:)
		!
		if( .not. velo_is_herm(A_ka))	then
			write(*,*)	'[gauge_trafo/DEBUG-MODE]: WARNING	(Hbar)-gauge A_ka IS NOT hermitian'	
		end if
		!
		if( .not. curv_is_herm(Om_kab) )	then
			write(*,*)	'[gauge_trafo/DEBUG-MODE]: WARNING	(Hbar)-gauge Om_kab IS NOT hermitian'				
		end if
		!
		!
		return
	end subroutine


	subroutine check_H_gauge_herm(kpt_rel, H_ka, A_ka, Om_kab, V_ka)
		real(dp),		intent(in)			::	kpt_rel(3)
		complex(dp),	intent(in)			::	H_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:),	&
												V_ka(:,:,:)
		character(len=30)					::	k_string
		!
		write(k_string,'(a,f6.2,a,f6.2,a,f6.2,a)')	'( ',kpt_rel(1),', ',kpt_rel(2),', ',kpt_rel(3),')'



		if(.not. velo_is_herm(A_ka)	)	then
			write(*,*)	"[get_wann_interp/DEBUG-MODE]:	WARNING (H)-gauge A_ka IS NOT hermitian at rel. kpt= "//k_string
		end if
		!
		!
		if(	.not. curv_is_herm( Om_kab)	) 	then
			write(*,*)	"[get_wann_interp/DEBUG-MODE]:	WARNING (H)-gauge Om_kab ARE NOT hermitian at rel. kpt= "//k_string
		end if
			!
		if(.not. velo_is_herm(H_ka)	)	then
			write(*,*)	"[get_wann_interp/DEBUG-MODE]:	WARNING (Hbar)-gauge H_ka IS NOT hermitian at rel. kpt= "//k_string
		end if
			!
		if(	velo_is_herm( V_ka )	) then
			write(*,*)	"[get_wann_interp/DEBUG-MODE]:	SUCCESS (H)-gauge V_ka IS hermitian at rel. kpt= "//k_string
		else
			write(*,*)	"[get_wann_interp/DEBUG-MODE]:	WARNING (H)-gauge V_KA IS NOT hermitian at rel. kpt= "//k_string
		end if
		!
		!
		return
	end subroutine







end module wann_interp