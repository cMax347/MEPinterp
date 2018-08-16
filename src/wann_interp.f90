module wann_interp
	!	module for Wannier interpolation
	!
	!	the interpolation scheme from 
	!			PRB 74, 195118 (2006) 
	!	was used
	use parameters,		only:		mpi_id,						&
									dp, i_dp, myExp,			&
									a_latt,	recip_latt, 		&
									my_Levi_Civita

	use matrix_math,	only:	zheevd_wrapper, matrix_comm






	implicit none


	private
	public					::		get_wann_interp	


	contains



!public:
	subroutine get_wann_interp(H_real, r_real, R_frac, kpt_rel, 	e_k, V_ka, A_ka, Om_ka  )
		complex(dp),					intent(in)				::	H_real(:,:,:)
		complex(dp),	allocatable, 	intent(inout)			::	r_real(:,:,:,:)
		real(dp),						intent(in)				::	R_frac(:,:), kpt_rel(3)	
		real(dp),						intent(out)				::	e_k(:)
		complex(dp),	allocatable,	intent(inout)			::	V_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout)			::	A_ka(:,:,:), Om_ka(:,:,:)
		
		complex(dp),	allocatable								::	U_k(:,:), H_ka(:,:,:), Om_kab(:,:,:,:), D_ka(:,:,:)

		!
		!
		allocate(	U_k(			size(H_real,1),size(H_real,2)	)			)		
		if(	allocated(V_ka)	)	then				
			allocate(	H_ka(	3,size(H_real,1),size(H_real,2)	 )		)
			!
			!
			if(	allocated(r_real)	) then
				allocate(	Om_kab(	3, 3, 	size(r_real,2), size(r_real,3) 	))
				allocate(	D_ka(		3,	size(r_real,2),	size(r_real,3)	))
			end if
		end if
		!
		!
		!ft onto k-space (W)-gauge
		call wann_interp_ft(H_real, r_real, R_frac, kpt_rel, U_k,  H_ka, A_ka, Om_kab)
		!get energies (H)-gauge
		call zheevd_wrapper(U_k, e_k)
		!
		!
		if(	allocated(V_ka)	)	then
			!do 	(W) -> (Hbar)
			call rotate_gauge(U_k, H_ka, 	A_ka, Om_kab )
			!
			!Velos (Hbar)->(H)
			call velo_gaugeTrafo(e_k, H_ka, A_ka, 	V_ka)


			!conn/curv	 (Hbar) -> (H)
			if( allocated(r_real) )	then
				call get_gauge_covar_deriv(e_k, H_ka, D_ka)
				!
				call conn_gaugeTrafo(D_ka, A_ka)
				call curv_gaugeTrafo(H_ka, A_ka, D_ka, Om_kab)
				!curv tens to vectors
				call om_tens_to_vect(Om_kab, Om_ka)
			end if
		end if
		!
		return
	end subroutine







!private
	subroutine wann_interp_ft(H_real, r_real, R_frac, kpt_rel, H_k,	H_ka, A_ka, Om_kab)			
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
		real(dp),						intent(in)				::	R_frac(:,:), kpt_rel(3)	
		complex(dp),					intent(out)				::	H_k(:,:)
		complex(dp),	allocatable,	intent(inout)			::	H_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp)												::	r_vect(3), kpt_abs(3)
		complex(dp)												::	ft_phase
		logical													::	use_pos_op, do_en_grad
		integer    												::	sc, a, b 
		!
		kpt_abs(1:3)	= 	matmul(		recip_latt(1:3,1:3)	, kpt_rel(1:3)	)	
		!
		H_k				=	dcmplx(0.0_dp)
		!
		!OPTIONAL energy gradients
		do_en_grad		= allocated(H_ka)
		if(do_en_grad)	H_ka	=	dcmplx(0.0_dp)
		!OPTIONAL position operator
		use_pos_op		= allocated(A_ka) .and. allocated(r_real) .and. allocated(Om_kab)
		if(use_pos_op)	then
			A_ka 		=	dcmplx(0.0_dp)
			Om_kab		=	dcmplx(0.0_dp)
		end if
		!			
		!
		!sum real space cells
		do sc = 1, size(R_frac,2)
			r_vect(:)	=	matmul(	a_latt(:,:),	R_frac(:,sc) )
			ft_phase	= 	myExp(	dot_product(kpt_abs(1:3),	r_vect(1:3)	)		)
			!
			!
			!Hamilton operator
			H_k(:,:)				= 	H_k(:,:)		+	ft_phase 					* H_real(:,:,sc)	
			!
			!
			do a = 1, 3
				!OPTIONAL energy gradients
				if( do_en_grad)		then
					H_ka(a,:,:) 		=	H_ka(a,:,:)		+	ft_phase * i_dp * r_vect(a) * H_real(:,:,sc)
				end if
				!OPTIONAL position operator
				if( use_pos_op )	then
					!connection
					A_ka(a,:,:)			=	A_ka(a,:,:)		+	ft_phase					* r_real(a,:,:,sc)
					!curvature
					do b = 1, 3
						Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:) + 	ft_phase * i_dp * r_vect(a) * r_real(b,:,:,sc)
						Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:) - 	ft_phase * i_dp * r_vect(b) * r_real(a,:,:,sc)
					end do
				end if 
			end do
			!
			!
		end do		
		!
		!
		return
	end subroutine




!gauge TRAFOs
	subroutine velo_gaugeTrafo(e_k, H_ka, A_ka, V_k)
		!	calc the (H)-gauge velocity matrix
		!
		!	PRB 74, 195118 (2006) EQ.(31)
		real(dp),						intent(in)		::		e_k(:)
		complex(dp),					intent(inout)	::		H_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout) 	::		A_ka(:,:,:)
		complex(dp),					intent(out)		::		V_k(:,:,:)
		integer											::		a, m, n
		!
		V_k	=	dcmplx(0.0_dp)
		!
		do a = 1, 3
			V_k(a,:,:)	=	V_k(a,:,:) +	H_ka(a,:,:)
			!
			if( allocated(A_ka)	) then
				do n = 1, size(V_k,3)
					do m = 1, size(V_k,2)
						V_k(a,m,n)	= V_k(a,m,n)	-	i_dp	*	(	e_k(m) - e_k(n)	)	*	A_ka(a,m,n)
					end do
				end do
			end if
			!
		end do
		!
		return
	end subroutine


	subroutine get_gauge_covar_deriv(e_k, H_ka,	D_ka )
		!	PRB 74, 195118 (2006)	EQ.(24)
		real(dp),			intent(in)		::	e_k(:)
		complex(dp),		intent(in)		::	H_ka(:,:,:)
		complex(dp),		intent(out)		::	D_ka(:,:,:)
		integer								::	m, n
		real(dp)							::	eDiff
		!
		D_ka(:,:,:)	=	dcmplx(0.0_dp)
		!
		do m = 1, size(D_ka,3)
			do n = 1, size(D_ka,2)
				if(	n/=	m )	then
					eDiff		=	e_k(m)	- e_k(n)
					if(abs(eDiff) < 1e-14_dp)	then
						eDiff	= sign(1e-14_dp,eDiff)
						write(*,'(a,i3,a)',advance="no")	'[#',mpi_id,';get_gauge_covar_deriv]:'
						write(*,'(a,i6,a,i6)')	' WARNING degenerate bands detetected n=',n,' m=',m
					end if
					!
					!
					D_ka(1:3,n,m)	=	H_ka(1:3,n,m) / 	eDiff
				end if
			end do
		end do 
		!
		return
	end subroutine




	subroutine conn_gaugeTrafo( D_ka, A_ka)
		!	PRB 74, 195118 (2006)	EQ.(25)
		complex(dp),		intent(in)		::	D_ka(:,:,:)
		complex(dp),		intent(inout)	::	A_ka(:,:,:)
		integer								::	a
		!
		do a = 1, 3
			A_ka(a,:,:)	=	A_ka(a,:,:)		+	i_dp	*	D_ka(a,:,:)
		end do
		!
		return
	end subroutine


	subroutine curv_gaugeTrafo(H_ka, A_ka, D_ka, Om_kab)
		!	PRB 74, 195118 (2006)	EQ.(27)
		complex(dp),		intent(in)		::	H_ka(:,:,:), A_ka(:,:,:), D_ka(:,:,:)
		complex(dp),		intent(inout)	::	Om_kab(:,:,:,:)
		complex(dp),	allocatable			::	mat_comm(:,:)
		integer								::	a, b
		!
		allocate(	mat_comm(	size(Om_kab,3),size(Om_kab,4)	)		)
		!
		do b = 1, 3
			do a = 1, 3
				call matrix_comm( D_ka(a,:,:), 	A_ka(b,:,:),		mat_comm(:,:)	)
				Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)		-			mat_comm(:,:)
				!
				!
				call matrix_comm( D_ka(b,:,:), 	A_ka(a,:,:),		mat_comm(:,:)	)
				Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)		+			mat_comm(:,:)
				!
				!
				call matrix_comm( D_ka(a,:,:), 	D_ka(b,:,:),		mat_comm(:,:)	)
				Om_kab(a,b,:,:)	=	Om_kab(a,b,:,:)		-	i_dp *	mat_comm(:,:)
			end do
		end do
		!
		!
		return 
	end subroutine

	



	subroutine rotate_gauge(U_k, H_ka, A_ka, Om_kab)
		complex(dp),					intent(in)		::	U_k(:,:)
		complex(dp), 					intent(inout)	::	H_ka(:,:,:)
		complex(dp), 	allocatable,	intent(inout)	::	A_ka(:,:,:), Om_kab(:,:,:,:)
		integer											::	a, b
		!
		do a = 1, 3
			call gauge_trafo(U_k,	H_ka(a,:,:))
			if( allocated(A_ka)		)	call gauge_trafo( U_k, A_ka(a,:,:))

			do b = 1,3 
				if( allocated(Om_kab)	)	call gauge_trafo( U_k, Om_kab(a,b,:,:))
			end do
		end do
		!
		return
	end subroutine


	subroutine gauge_trafo(U_mat, M_mat)
		!	does a gauge trafo of the gauge-covariant matrix M_mat
		!	by applying the unitary matrix U_mat as shown in 
		!	Eq.(21), PRB 74, 195118, (2006)
		!
		!	since U is defined differently in the paper and in the lapack
		!	consider the following
		!
		!	want: 	M^(hbar) 	
		!	
		!	 paper:		 U^*	H^(W)	U 	= H^(H)
		!			and	M^(hbar) = U^* M^(W) U
		!
		!	lapack:
		!			 U	H^(W)	U^* 	= H^(H)	
		!	therefore here M^(hbar) = U M^(W) U^*

		complex(dp),		intent(in)		::	U_mat(:,:)
		complex(dp),		intent(inout)	::	M_mat(:,:)
		complex(dp),	allocatable			::	U_dag(:,:)
		!
		allocate(		U_dag(	size(U_mat,1),size(U_mat,2)	)			)
		!
		U_dag 	= conjg(	transpose( U_mat )	)
		!
		M_mat	= matmul(	M_mat	,	U_dag	)
		M_mat	= matmul(	U_mat	,	M_mat	)
		!
		return
	end subroutine












!helpers:
	subroutine om_tens_to_vect(om_tens, om_vect)
		!	converts om_tens to a vector by applying Levi Cevita Tensor
		!	see PRB 74, 195118 (2006)	Eq.(5)
		!
		complex(dp),		intent(in)		::	om_tens(:,:,:,:)
		complex(dp),		intent(out)		::	om_vect(:,:,:)
		integer								::	a, b, c
		!
		om_vect	=	dcmplx(0.0_dp)
		!
		do c = 1, 3
			!
			do b = 1, 3
				do a = 1,3
					om_vect(c,:,:)	= om_vect(c,:,:) + real(my_Levi_Civita(a,b,c),dp) * om_tens(a,b,:,:)
				end do
			end do
			!
		end do
		return

	end subroutine





end module wann_interp