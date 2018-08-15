module wann_interp
	use parameters,		only:		mpi_id,						&
									dp, i_dp, myExp,			&
									a_latt,	recip_latt, 		&
									do_gauge_trafo

	use matrix_math,	only:	zheevd_wrapper






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
		complex(dp),	allocatable,	intent(out)				::	V_ka(:,:,:)
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
			call rotate_gauge(U_k, 	A_ka, Om_kab, H_ka )
			!
			!Velos (Hbar)->(H)
			call velo_gaugeTrafo(e_k, H_ka, A_ka, 	V_ka)



			!conn/curv	 (Hbar) -> (H)
			if( allocated(r_real) )	then
				call conn_gaugeTrafo(H_ka, A_ka, D_ka)
				call curv_gaugeTrafo(H_ka, A_ka, D_ka, Om_kab)
				!curv tens to vectors
				call om_tens_to_vect(Om_kab, Om_ka)
			end if
		end if


		return
	end subroutine







!private
	subroutine velo_gaugeTrafo(e_k, H_ka, A_ka, V_k)
		real(dp),						intent(in)		::		e_k(:)
		complex(dp),					intent(inout)	::		H_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout) 	::		A_ka(:,:,:)
		complex(dp),					intent(out)		::		V_k(:,:,:)
		integer											::		x, m, n
		!
		V_k	=	dcmplx(0.0_dp)
		!
		!
		do x = 1, 3
			V_k(x,:,:)	=	V_k(x,:,:) +	H_ka(x,:,:)
			!
			if( allocated(A_ka)	) then
				do n = 1, size(V_k,3)
					do m = 1, size(V_k,2)
						V_k(x,m,n)	= V_k(x,m,n)	-	i_dp	*	(	e_k(m) - e_k(n)	)	*	A_ka(x,m,n)
					end do
				end do
			end if
			!
		end do
		!
		!
		return
	end subroutine


	subroutine conn_gaugeTrafo(H_ka, A_ka, D_ka)
		complex(dp),		intent(in)		::	H_ka(:,:,:)
		complex(dp),		intent(inout)	::	A_ka(:,:,:)
		complex(dp),		intent(out)		::	D_ka(:,:,:)

		D_ka	= 	dcmplx(0.0_dp)

		write(*,*)	"WARNING [wann_interp/do_conn_gaugeTrafo] is not implemented yet"
		return
	end subroutine

	subroutine curv_gaugeTrafo(H_ka, A_ka, D_ka, Om_kab)
		complex(dp),		intent(in)		::	H_ka(:,:,:), A_ka(:,:,:), D_ka(:,:,:)
		complex(dp),		intent(inout)	::	Om_kab(:,:,:,:)

		write(*,*)	"WARNING [wann_interp/do_curv_gaugeTrafo] is not implemented yet"
		return 
	end subroutine

	


	subroutine rotate_gauge(U_k, A_ka, Om_kab, H_ka)
		complex(dp),					intent(in)		::	U_k(:,:)
		complex(dp), 					intent(inout)	::	H_ka(:,:,:)
		complex(dp), 	allocatable,	intent(inout)	::	A_ka(:,:,:), Om_kab(:,:,:,:)
		integer											::	a, b
		
		write(*,*)	"[wann_interp/rotate_gauge]:	WARNING not implemented yet"
		!
		if( allocated(A_ka)		)	then

		end if
		!
		if( allocated(Om_kab)	)	then

		end if
		!
		return
	end subroutine



	subroutine wann_interp_ft(H_real, r_real, R_frac, kpt_rel, H_k,	H_ka, A_ka, Om_kab)			
		!	interpolates real space Ham and position matrix to k-space,
		!	according to
		!		PRB 74, 195118 (2006)
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
			do a = 1, 3
				do b = 1,3
					om_vect(c,:,:)	= om_vect(c,:,:) + real(levi_cicvita(a,b,c),dp) * om_tens(a,b,:,:)
				end do
			end do
			!
		end do
		return

	end subroutine

	integer function levi_cicvita(a,b,c)
		integer							::	a, b, c
		!
		levi_cicvita = 0
		!
		if(			(a==1 .and. b==2 .and. c==3) .or. (a==2 .and. b==3 .and. c==1) .or. (a==3 .and. b==1 .and. c==2)	) then
			levi_cicvita = 1
		else if(	(a==3 .and. b==2 .and. c==1) .or. (a==1 .and. b==3 .and. c==2) .or.	(a==2 .and. b==1 .and. c==3)	) then
			levi_cicvita = -1
		end if
		!
		return
	end function



	subroutine gauge_trafo(U_mat, M_mat)
		!	does a gauge trafo of the gauge-covariant matrix M_mat
		!	by applying the unitary matrix U_mat as shown in 
		!	Eq.(21), PRB 74, 195118, (2006)
		complex(dp),		intent(in)		::	U_mat(:,:)
		complex(dp),		intent(inout)	::	M_mat(:,:)
		complex(dp),	allocatable			::	U_dag(:,:)
		!
		allocate(		U_dag(	size(U_mat,1),size(U_mat,2)	)			)
		!
		U_dag = 	conjg(		transpose( U_mat )		)
		!
		M_mat	= matmul(	M_mat	,	U_mat	)
		M_mat	= matmul(	U_dag	,	M_mat	)
		!
		return
	end subroutine





end module wann_interp