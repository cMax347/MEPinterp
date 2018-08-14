module wann_interp
	use parameters,		only:		mpi_id,						&
									dp, i_dp, myExp,			&
									a_latt,	recip_latt, 		&
									do_gauge_trafo








	implicit none


	private
	public					::		wann_interp_ft,		velo_interp		


	contains



!public:
	subroutine wann_interp_ft(H_real, r_real, R_frac, kpt_rel, H_k,	H_ka, A_ka, Om_ka)			
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
		complex(dp),	allocatable,	intent(inout)			::	H_ka(:,:,:), A_ka(:,:,:), Om_ka(:,:,:)
		complex(dp),	allocatable								::	Om_kab(:,:,:,:)
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
		use_pos_op		= allocated(A_ka) .and. allocated(r_real) .and. allocated(Om_ka)
		if(use_pos_op)	then
			allocate(	Om_kab(	3, 3, size(r_real,2), size(r_real,3) ))
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
		!get curvature in desired vector convention
		if(use_pos_op)	call om_tens_to_vect(Om_kab, Om_ka)
		!
		return
	end subroutine



	subroutine velo_interp(U_k, en_k, H_ka, A_ka, V_k)
		complex(dp),					intent(in)		::		U_k(:,:)
		real(dp),						intent(in)		::		en_k(:)
		complex(dp),					intent(inout)	::		H_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout) 	::		A_ka(:,:,:)
		complex(dp),					intent(out)		::		V_k(:,:,:)
		integer											::		x, m, n
		!
		V_k	=	dcmplx(0.0_dp)
		!
		!
		do x = 1, 3

			if( do_gauge_trafo)		call gauge_trafo(U_k,	H_ka(x,:,:)		)

			V_k(x,:,:)	=	V_k(x,:,:) +	H_ka(x,:,:)

			if( allocated(A_ka)	) then
				if( do_gauge_trafo)		call gauge_trafo(U_k,	A_ka(x,:,:) )
				!
				do n = 1, size(V_k,3)
					do m = 1, size(V_k,2)
						V_k(x,m,n)	= V_k(x,m,n)	-	i_dp	*	(	en_k(m) - en_k(n)	)	*	A_ka(x,m,n)
					end do
				end do
			end if
			!
			!
		end do


		return
	end subroutine









!private:
	subroutine om_tens_to_vect(om_tens, om_vect)
		complex(dp),		intent(in)		::	om_tens(:,:,:,:)
		complex(dp),		intent(out)		::	om_vect(:,:,:)

		om_vect	=	dcmplx(0.0_dp)
		write(*,*)	'[om_tens_to_vect]: WARNING NOT IMPLEMENTED YET, SET OM_VECT TO ZERO'
		return
	end subroutine



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