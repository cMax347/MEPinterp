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
									blas_matmul,			&
									is_herm_mat,			&
									is_skew_herm_mat		
	use input_paras,	only:		kubo_tol, debug_mode
	use wrapper_3q,		only:		get_ham






	implicit none


	private
	public					::		get_wann_interp, W_to_H_gaugeTRAFO	


	contains









!public:
	subroutine get_wann_interp(		do_gauge_trafo,				&
									kpt_idx, kpt_rel, 			&
									e_k, V_ka, 					&
									dummy_conn, dummy_curv		&
							)
		!
		!	interpolates the k-space:
		!			-	U_k	: 		hamiltonian, zheevd_wrapper overwrites to unitary gauge_trafo matrix ("hamiltonian gauge")
		!			-	dummy_conn:		since setup is in k-space no connection is given	
		!			-	dummy_curv:		since setup is in k-space no curv is given
		!
		!
		!	see 	>>>	PRB 74, 195118 (2006)	<<<			for more details on Wannier interpolation
		!
		logical,						intent(in)				::	do_gauge_trafo
		integer,						intent(in)				::	kpt_idx
		real(dp),						intent(in)				::	kpt_rel(3)															
		real(dp),		allocatable,	intent(inout)			::	e_k(:)
		complex(dp),	allocatable,	intent(inout)			::	V_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout)			::	dummy_conn(:,:,:), dummy_curv(:,:,:,:)	!DO NOT ALLOCATE DUMMYS, DUMMY
		complex(dp),	allocatable								::	U_k(:,:)
		!
		!


		!todo:	get k_space hamiltonian and store it in U_k
		write(*,*)	"[get_wann_interp]	WARNING	get_ham takes relative kpt for now (todo: fix this)"	
		call get_ham(kpt_rel,	U_k, V_ka)



		!get energies (H)-gauge
		if(.not. allocated(e_k))	allocate(e_k(size(U_k,1)))
		call zheevd_wrapper(U_k, e_k)
		!
		!rotate back to (H)-gauge
		if( allocated(V_ka)	.and. do_gauge_trafo	)			call W_to_H_gaugeTRAFO(e_k, U_k, V_ka, dummy_conn, dummy_curv)
		!
		!
		!	DEBUG
		if(debug_mode)	call check_H_gauge_herm(kpt_idx, kpt_rel, dummy_conn, dummy_curv, V_ka)
		!
		return
	end subroutine







!private:
!	subroutine FT_R_to_k(H_real, r_real, a_latt, recip_latt, R_frac, kpt_rel, H_k,	H_ka, A_ka, Om_kab)			
!		!	interpolates real space Ham and position matrix to k-space,
!		!	according to
!		!		PRB 74, 195118 (2006)		EQ.(37)-(40)
!		!
!		!
!		!	->	only the H_real, and H_k have to be allocated
!		!	->	all other quantities are only calculated if allocated
!		!
!		complex(dp),					intent(in)				::	H_real(:,:,:)
!		complex(dp),	allocatable, 	intent(inout)			::	r_real(:,:,:,:)
!		real(dp),						intent(in)				::	a_latt(3,3), recip_latt(3,3),	&
!																	R_frac(:,:), kpt_rel(3)	
!		complex(dp),					intent(out)				::	H_k(:,:)
!		complex(dp),	allocatable,	intent(inout)			::	H_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:)
!		real(dp)												::	R_abs(3), kpt_abs(3), ft_angle
!		complex(dp)												::	ft_phase
!		logical													::	use_pos_op, do_en_grad
!		integer    												::	sc, a, b 
!		!
!		!jobs
!		do_en_grad		= allocated(H_ka)
!		use_pos_op		= allocated(A_ka) .and. allocated(r_real) .and. allocated(Om_kab)
!		!
!		!get absolute kpt
!		kpt_abs		= 	matmul(	recip_latt	, kpt_rel	)	
!		!
!		!init
!						H_k		= 0.0_dp
!		if(do_en_grad)	H_ka	= 0.0_dp
!		if(use_pos_op)	A_ka	= 0.0_dp
!		if(use_pos_op)	Om_kab	= 0.0_dp	
!		!
!		!
!		!sum real space cells
!		do sc = 1, size(R_frac,2)
!			R_abs(:)	=	matmul(	a_latt(:,:),	R_frac(:,sc) )
!			ft_angle	=	dot_product(kpt_abs(1:3),	R_abs(1:3))
!			ft_phase	= 	cmplx(	cos(ft_angle), sin(ft_angle)	,	dp	)
!			!
!			!ft_phase	=	cmplx(1.0_dp, 0.0_dp, dp)
!			!
!			!Hamilton operator
!			H_k(:,:)				= 	H_k(:,:)			+	ft_phase 					* H_real(:,:,sc)	
!			!
!			do a = 1, 3
!				!OPTIONAL energy gradients
!				if( do_en_grad)		then
!					H_ka(a,:,:) 		=	H_ka(a,:,:)		+	ft_phase * i_dp * R_abs(a) * H_real(:,:,sc)
!				end if
!				!OPTIONAL position operator
!				if( use_pos_op )	then
!					!connection
!					A_ka(a,:,:)			=	A_ka(a,:,:)		+	ft_phase					* r_real(a,:,:,sc)
!					!curvature
!					do b = 1, 3
!						Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:) + 	ft_phase * i_dp * R_abs(a) * r_real(b,:,:,sc)
!						Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:) - 	ft_phase * i_dp * R_abs(b) * r_real(a,:,:,sc)
!					end do
!				end if 
!			end do
!			!
!			!
!		end do		
!		!
!		!
!		if(debug_mode) 	then
!			call check_ft_phase(a_latt, R_frac, kpt_abs)
!			call check_W_gauge_herm(kpt_rel, H_k, H_Ka, A_ka, Om_kab)
!		end if
!		!
!		return
!	end subroutine



	pure subroutine get_velo(e_k, H_ka, A_ka, V_ka)
		!	calc the (H)-gauge velocity matrix
		!
		!	PRB 74, 195118 (2006) EQ.(31)
		real(dp),						intent(in)		::		e_k(:)
		complex(dp),					intent(in)		::		H_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout) 	::		A_ka(:,:,:)
		complex(dp),					intent(out)		::		V_ka(:,:,:)
		complex(dp)										::		eDiff
		integer											::		m, n
		!
		V_ka		=	H_ka
		!
		!
		if( allocated(A_ka)	) then
			do m = 1, size(V_ka,3)
				do n = 1, size(V_ka,2)
					if(	n >	m )	then
						eDiff	=	cmplx(		e_k(m) - e_k(n),		0.0_dp,	dp)
						!
						V_ka(:,n,m)	= V_ka(:,n,m)	- i_dp	* 	eDiff	* 	A_ka(:,n,m)
						V_ka(:,m,n)	= V_ka(:,m,n)	+ i_dp	*	eDiff	*	A_ka(:,m,n)
					end if
				end do
			end do
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
		!conn/curv	 (Hbar) -> (H)
		if( allocated(A_ka) )	then
			stop "W_to_H_gaugeTRAFO: in 3qstate model to connection is not used (and should not be allocated) "
			if(debug_mode)	call check_Hbar_gauge_herm(H_ka, A_ka, Om_kab)			! if(.not.allocated(A_ka)) then  (Hbar)-gauge is same as (H)-gauge (therefore avoid testing same twice)
			!
			!
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
		complex(dp),	allocatable						::	U_dag(:,:)
		integer											::	a, b
		!
		allocate(	U_dag(size(U_k,1),size(U_k,2))	)
		U_dag	=	conjg(	transpose(	U_k(:,:)	))

		do a = 1, 3
										H_ka(a,:,:)	=	blas_matmul(	blas_matmul(U_dag,	H_ka(a,:,:))	,	U_k	)
			if( allocated(A_ka)		)	A_ka(a,:,:)	=	blas_matmul(	blas_matmul(U_dag,	A_ka(a,:,:))	,	U_k	)
			if( allocated(Om_kab)	)then
				do b = 1,3 
									Om_kab(a,b,:,:)	=	blas_matmul(	blas_matmul(U_dag,	Om_kab(a,b,:,:)),	U_k	)
				end do
			end if
		end do
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
						write(*,'(a,i6,a,i6)',advance="no")		'WARNING degenerate bands detetected n=',n,' m=',m
						write(*,*)	"	(this might lead to precission issues with gauge covariant derivate D_ka)"
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
	subroutine	check_ft_phase(a_latt, R_frac, kpt_abs)
		real(dp),				intent(in)		::	a_latt(3,3), R_frac(:,:), kpt_abs(3)
		real(dp)								::	R_abs(3), ft_angle
		complex(dp)								::	ft_phase	
		integer									::	sc	
		!
		ft_phase		=	cmplx(0.0_dp, 0.0_dp, dp)
		do sc = 1, size(R_frac,2)
			R_abs(:)	=	matmul(	a_latt(:,:),	R_frac(:,sc) )
			ft_angle	=	dot_product(kpt_abs(1:3),	R_abs(1:3))
			ft_phase	=	ft_phase	+	cmplx(	cos(ft_angle), sin(ft_angle)	,	dp	)
		end do
		if(	aimag(ft_phase)	> 1e-6_dp	) write(*,*)	"[FT_R_to_k/DEBUG-MODE]: WARNING sum[ imag(ft_phase) ] =",aimag(ft_phase)," /=0"
		ft_phase		=	cmplx(0.0_dp, 0.0_dp, dp)	
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


	subroutine check_H_gauge_herm(kpt_idx, kpt_rel, A_ka, Om_kab, V_ka)
		integer,						intent(in)			::	kpt_idx
		real(dp),						intent(in)			::	kpt_rel(3)
		complex(dp),	allocatable,	intent(in)			::	A_ka(:,:,:), Om_kab(:,:,:,:), 	V_ka(:,:,:)
		real(dp)											::	max_err
		character(len=31)									::	k_string
		logical												::	conn, curv, velo
		!
		write(k_string,'(a,f6.2,a,f6.2,a,f6.2,a)')	'( ',kpt_rel(1),', ',kpt_rel(2),', ',kpt_rel(3),') '
		!
		!
		!	CONNECTION
		if(allocated(A_ka)) then
			conn	=	velo_is_herm(A_ka, max_err)
			if(.not. conn) 		write(*,'(a,f16.7)')	"[check_H_gauge_herm/DEBUG-MODE]:	"								//	&
															"WARNING (H)-gauge A_ka IS NOT hermitian at rel. kpt= "			//	&
															k_string//"max_err=", max_err
		else
			conn	=	.true.
		end if
		!
		!
		!	CURVATURE
		if(allocated(Om_kab)) then
			curv	=	curv_is_herm( Om_kab, max_err)
			if(.not. curv) 		write(*,'(a,f16.7)')	"[check_H_gauge_herm/DEBUG-MODE]:	"						 		//	&
															"WARNING (H)-gauge Om_kab IS NOT hermitian at rel. kpt= "		//	&
															k_string//"max_err=", max_err
		else
			curv	=	.true.
		end if
		!
		!
		!	VELOCITY
		if(allocated(V_ka)) then
			velo	=	velo_is_herm( V_ka, max_err )
			if(.not. velo) 		write(*,'(a,f16.7)')	"[check_H_gauge_herm/DEBUG-MODE]:	"								//	&
															"WARNING (H)-gauge V_KA IS NOT hermitian at rel. kpt= "			//	&
															k_string//"max_err=", max_err
		else
			velo	=	.true.
		end if
		!
		!
		!	SUCCESS MESSAGE
		if(	conn .and. curv .and. velo ) 	write(*,'(a,i8)')	"[check_H_gauge_herm/DEBUG-MODE]:	"						//	&
															"SUCCESS (H)-gauge quantities (conn,curv,velo) are hermitian "	//  &
															" at  #kpt= ", kpt_idx		
		!
		return
	end subroutine










end module wann_interp