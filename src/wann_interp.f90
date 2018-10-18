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
	use file_IO,		only:		write_velo





	implicit none


	private
	public					::		get_wann_interp	


	contains



!public:
	subroutine get_wann_interp(		do_gauge_trafo, H_real, r_real, 					&
									a_latt, recip_latt, R_frac, kpt_idx, kpt_rel, 		&
									e_k, V_ka, A_ka, Om_kab								&
							)
		logical,						intent(in)				::	do_gauge_trafo
		complex(dp),					intent(in)				::	H_real(:,:,:)
		complex(dp),	allocatable, 	intent(inout)			::	r_real(:,:,:,:)
		integer,						intent(in)				::	kpt_idx
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
		if(	allocated(V_ka)							)			call get_velo(e_k, H_ka, A_ka, 	V_ka)
		!
		!
		!	DEBUG
		if(debug_mode)	then
			call check_H_gauge_herm(kpt_rel, H_ka, A_ka, Om_kab, V_ka)
			if(allocated(V_ka))	call debug_write_velo_files(kpt_idx, e_k, H_ka, A_ka, V_ka)
		end if
		!
		return
	end subroutine


	subroutine debug_write_velo_files(kpt_idx, e_k, H_ka, A_ka, V_ka)
		integer,						intent(in)		::	kpt_idx
		real(dp),						intent(in)		::	e_k(:)
		complex(dp),					intent(in)		::	H_ka(:,:,:), V_ka(:,:,:)
		complex(dp), 	allocatable, 	intent(in)		::	A_ka(:,:,:)
		complex(dp),	allocatable						::	Vtemp(:,:,:)
		character(len=10)								::	fname
		character(len=120)								::	info_string
		integer											::	n, m
		!
		!
		if(allocated(A_ka))	 then
			write(fname,*)		"velo_Aka"
			write(info_string,*)	"#	Berry Connection contribution to velocities: V_ka(n,m) =  - i (E_m - E_n) A_ka(n,m)"
			!
			allocate(		Vtemp(3, size(A_ka,2), size(A_ka,3))			)
			!
			Vtemp	= cmplx(0.0_dp, 0.0_dp,dp)
			do n = 1, size(A_ka,2)
				do m = 1, size(A_ka,3)
					Vtemp(:,n,m)	=	-i_dp	* cmplx( e_k(m)	- e_k(n), 0.0_dp, dp) *	A_ka(:,n,m)
				end do
			end do
			call write_velo(kpt_idx, fname, info_string, Vtemp)
			!
			!
			write(fname,*)		"velo_Vka"
			write(info_string,*)	"#	(H)-gauge velo velocity: V_ka(n,m) =  H_ka(n,m)  - i (E_m - E_n) A_ka(n,m)"
			call write_velo(kpt_idx, fname, info_string, V_ka)
		end if
		!
		write(fname,*)		"velo_Hka"
		write(info_string,*)	"#	Hamiltonian contribution to velocities: V_ka =  H_ka"
		call write_velo(kpt_idx, fname, info_string, H_ka)
		!
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
		!	test 
		ft_phase	=	cmplx(0.0_dp, 0.0_dp, dp)
		do sc = 1, size(R_frac,2)
			R_abs(:)	=	matmul(	a_latt(:,:),	R_frac(:,sc) )
			ft_angle	=	dot_product(kpt_abs(1:3),	R_abs(1:3))
			ft_phase	=	ft_phase	+	cmplx(	cos(ft_angle), sin(ft_angle)	,	dp	)
		end do
		if(	aimag(ft_phase)	> 1e-6_dp	) write(*,*)	"[FT_R_to_k]: Warning sum[ imag(ft_phase) ] =",aimag(ft_phase)," /=0"

		!
		!sum real space cells
		do sc = 1, size(R_frac,2)
			R_abs(:)	=	matmul(	a_latt(:,:),	R_frac(:,sc) )
			ft_angle	=	dot_product(kpt_abs(1:3),	R_abs(1:3))
			ft_phase	= 	cmplx(	cos(ft_angle), sin(ft_angle)	,	dp	)
			!
			!ft_phase	=	cmplx(1.0_dp, 0.0_dp, dp)
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
				if( .not. is_skew_herm_mat(D_ka(a,:,:), max_err)	)	then
					write(*,'(a,i1,a,f16.7)')	"[get_gauge_covar_deriv/DEBUG-MODE]:	WARNING D_(k,a=",a,&
																						") is not skew hermitian, max_err=", max_err
				
					!if( .not. is_herm_mat(D_ka(a,:,:),max_err))	then
					!	write(*,'(a,i1,a,f16.7)')	"[get_gauge_covar_deriv/DEBUG-MODE]:	WARNING D_(k,a=",a,&
					!																	") is neither hermitian, max_err=", max_err
					!end if
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



!DEBUG helpers
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
		logical								::	basis, velo
		!
		basis	=	size(V_ka,2) == size(V_ka,3) 	
		!
		!
		if(basis)	then
			velo_is_herm	=	.true.
			do a = 1, 3
				if(.not. is_herm_mat(V_ka(a,:,:),max_err)) 	then
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
		real(dp),						intent(in)		::		kpt_rel(3)
		complex(dp),					intent(in)		::		H_k(:,:)
		complex(dp),	allocatable, 	intent(in)		::		H_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:) 
		real(dp)										::		max_err	
		character(len=30)								::		k_string
		!
		write(k_string,'(a,f6.2,a,f6.2,a,f6.2,a)')	'( ',kpt_rel(1),', ',kpt_rel(2),', ',kpt_rel(3),')'
		!
		!
		if(	.not. is_herm_mat(H_k,	max_err) ) then
			write(*,'(a,f16.7,a,200(f6.2))')	"[check_w_gauge_herm/DEBUG-MODE]:	WARNING (W)-gauge H_k IS NOT hermitian(max_err=",&
																				max_err,") at rel. kpt="//k_string
		end if
		!
		if(allocated(H_ka)) then
			if( .not. velo_is_herm(H_ka,	max_err)) then
				write(*,'(a,f16.7,a,200(f6.2))')	"[check_w_gauge_herm/DEBUG-MODE]:	WARNING (W)-gauge H_ka IS NOT hermitian(max_err=",&
																					max_err,") at rel. kpt="//k_string
			end if
		end if
		!
		!
		if(allocated(A_ka)) then
			if( .not. velo_is_herm( A_ka,	max_err)	) 	then
				write(*,*)						"[check_w_gauge_herm/DEBUG-MODE]:	WARNING (W)-gauge A_ka IS NOT hermitian(max_err=",&
																					max_err,") at rel. kpt= "//k_string
			end if
		end if
		!
		if(allocated(Om_kab))	then
			if( .not. curv_is_herm( Om_kab,	max_err)	)	then
				write(*,*)						"[check_w_gauge_herm/DEBUG-MODE]:	WARNING (W)-gauge OM_kab IS NOT hermitian(max_err=",&
																					max_err,") at rel. kpt= "//k_string
			end if
		end if
		!
		!
		return
	end subroutine


	subroutine check_Hbar_gauge_herm(A_ka, Om_kab)
		complex(dp), allocatable,	intent(in)			::	A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp)							::	max_err
		!
		if(allocated(A_ka)) then
			if( .not. velo_is_herm(A_ka,max_err))	then
				write(*,'(a,f16.7)')	'[gauge_trafo/DEBUG-MODE]: WARNING	(Hbar)-gauge A_ka IS NOT hermitian, max_err=',max_err	
			end if
		end if
		!
		if(allocated(Om_kab))	then
			if( .not. curv_is_herm(Om_kab,max_err) )	then
				write(*,'(a,f16.7)')	'[gauge_trafo/DEBUG-MODE]: WARNING	(Hbar)-gauge Om_kab IS NOT hermitian, max_err=',max_err				
			end if
		end if
		!
		!
		return
	end subroutine


	subroutine check_H_gauge_herm(kpt_rel, H_ka, A_ka, Om_kab, V_ka)
		real(dp),						intent(in)			::	kpt_rel(3)
		complex(dp),	allocatable,	intent(in)			::	H_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:), 	V_ka(:,:,:)
		real(dp)											::	max_err
		character(len=30)									::	k_string
		!
		write(k_string,'(a,f6.2,a,f6.2,a,f6.2,a)')	'( ',kpt_rel(1),', ',kpt_rel(2),', ',kpt_rel(3),')'


		if(allocated(A_ka)) then
			if(.not. velo_is_herm(A_ka, max_err)	)	then
				write(*,'(a,f16.7)')	"[get_wann_interp/DEBUG-MODE]:	WARNING (H)-gauge A_ka IS NOT hermitian at rel. kpt= "//	&
																								k_string//"max_err=", max_err
			end if
		end if
		!
		!
		if(allocated(Om_kab)) then
			if(	.not. curv_is_herm( Om_kab, max_err)	) 	then
				write(*,'(a,f16.7)')	"[get_wann_interp/DEBUG-MODE]:	WARNING (H)-gauge Om_kab ARE NOT hermitian at rel. kpt= "//		&
																								k_string//"max_err=", max_err
			end if
		end if
			!
		if(allocated(H_ka)) then
			if(.not. velo_is_herm(H_ka, max_err)	)	then
				write(*,'(a,f16.7)')	"[get_wann_interp/DEBUG-MODE]:	WARNING (Hbar)-gauge H_ka IS NOT hermitian at rel. kpt= "//		&
																								k_string//"max_err=", max_err
			end if
		end if
			!
		if(allocated(V_ka)) then
			if(	velo_is_herm( V_ka, max_err )	) then
				write(*,*)	"[get_wann_interp/DEBUG-MODE]:	SUCCESS (H)-gauge V_ka IS hermitian at rel. kpt= "//k_string
			else
				write(*,'(a,f16.7)')	"[get_wann_interp/DEBUG-MODE]:	WARNING (H)-gauge V_KA IS NOT hermitian at rel. kpt= "//		&
																								k_string//"max_err=", max_err
			end if
		end if
		!
		!
		return
	end subroutine







end module wann_interp