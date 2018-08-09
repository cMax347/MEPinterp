module wann_interp
	use parameters,		only:		mpi_id,						&
									dp, i_dp, myExp,			&
									a_latt,	recip_latt, 		&
									do_gauge_trafo


	implicit none


	private
	public					::		wann_interp_ft,		velo_interp		


	contains



	subroutine wann_interp_ft(H_real, r_real, R_frac, kpt_rel, H_k,	H_ka, A_ka	)										!H_real, r_mat,  R_vect,	kpt_rel_latt,	U_k,	H_ka, A_ka)	
		!	interpolates real space Ham and position matrix to k-space,
		!	according to
		!		PRB 74, 195118 (2006)
		!
		complex(dp),					intent(in)				::	H_real(:,:,:)
		complex(dp),	allocatable, 	intent(inout)			::	r_real(:,:,:,:)
		real(dp),						intent(in)				::	R_frac(:,:), kpt_rel(3)	
		complex(dp),					intent(out)				::	H_k(:,:), H_ka(:,:,:)
		complex(dp),	allocatable,	intent(inout)			::	A_ka(:,:,:)
		real(dp)												::	r_vect(3), kpt_abs(3)
		complex(dp)												::	ft_phase
		logical													::	do_pos
		integer    												::	sc
		!
		do_pos	= allocated(	A_ka	) .and. allocated(	r_real	)
		!
		
		!
		H_k				=	dcmplx(0.0_dp)
		H_ka			=	dcmplx(0.0_dp)
		if(do_pos)	A_ka=	dcmplx(0.0_dp)
		!
		!
		kpt_abs(1:3)	= 	matmul(		, kpt_rel(1:3)	)					!todo


		!
		!SUM OVER REAL SPACE CELLS
		do sc = 1, size(R_frac,2)
			r_vect(:)	=	matmul(	a_latt(:,:),	R_frac(:,sc) )
			ft_phase	= 	myExp(	dot_product(kpt_abs(1:3),	r_vect(1:3)	)		)

!			write(*,*)	"[#",mpi_id,"; wann_interp_ft]:	now start with ft at kpt_rel=(",	(kpt_rel(x),", ", x=1,3   )		,	") and R_vect=(",(r_vect(x),", ", x=1,3   ),")."	
			write(*,'(a,i3,a)',			advance="no")	"[#",mpi_id,"; wann_interp_ft]:	ft_phase=("
			write(*,'(f6.2,a,f6.2)',	advance="no")	dreal(ft_phase),'+ i ',dimag(ft_phase)
			write(*,'(a,i3,a)',			advance="no")	") sc=",sc," kpt_rel=("
			write(*,'(f6.2,a,f6.2,a,f6.2,a)')			kpt_rel(1),", ",kpt_rel(2),", ",kpt_rel(3),")."

			!hamilton operator
			H_k(:,:)		= 	H_k(:,:)		+	ft_phase 					* H_real(:,:,sc)
			H_ka(1,:,:)		= 	H_ka(1,:,:)		+	ft_phase * i_dp * r_vect(1)	* H_real(:,:,sc)	
			H_ka(2,:,:)		= 	H_ka(2,:,:)		+	ft_phase * i_dp * r_vect(2)	* H_real(:,:,sc)	
			H_ka(3,:,:)		= 	H_ka(3,:,:)		+	ft_phase * i_dp * r_vect(3)	* H_real(:,:,sc)	

			!positional operator
			if( do_pos ) then
				A_ka(1,:,:)	=	A_ka(1,:,:)		+	ft_phase 					* r_real(1,:,:,sc)	
				A_ka(2,:,:)	=	A_ka(2,:,:)		+	ft_phase 					* r_real(2,:,:,sc)	
				A_ka(3,:,:)	=	A_ka(3,:,:)		+	ft_phase 					* r_real(3,:,:,sc)	
			end if
		end do		



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