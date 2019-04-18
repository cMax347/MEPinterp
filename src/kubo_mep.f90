module kubo_mep
	!
	use constants,		only:	dp, i_dp, pi_dp
	use matrix_math,	only:	get_levi_civita,		& 
								convert_tens_to_vect
	use statistics,		only:	fd_stat
	implicit none
	!
	private
	public			::			kubo_mep_CS,		&
								kubo_mep_LC,		&
								kubo_mep_IC
							
contains

	pure function kubo_mep_CS(A_ka, Om_kab, fd_distrib) result(cs_tens)
		complex(dp),	allocatable, 	intent(in)		::	A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp),						intent(in)		::	fd_distrib(:,:)
		real(dp),		allocatable						::	cs_tens(:,:,:), fd_lst(:)
		complex(dp),	allocatable						::	cs_scal(:)
		complex(dp)										::	om_vect(3)
		integer											::	n0, i, n_ef
		!
		n_ef	=	size(fd_distrib,1)
		allocate(cs_tens(3,3,n_ef))
		cs_tens	=	0.0_dp
		!
		if(	allocated(A_ka)		.and.		allocated(Om_kab)		)then
			allocate(	fd_lst(		n_ef	))
			allocate(	cs_scal(	n_ef	))
			cs_scal	= 0.0_dp
			!
			!	sum over valence
			do n0 = 1, size(Om_kab,3)
				call convert_tens_to_vect(	Om_kab(:,:,n0,n0),		om_vect(:)	)
				!
				cs_scal(:)	= 	cs_scal(:) 	+	 0.5_dp  * 	fd_distrib(:,n0) 							&
												* dreal(	dot_product( A_ka(:,n0,n0)	, om_vect(:) ))
			end do
			!
			do i = 1, 3
				cs_tens(i,i,:)	= real(cs_scal(:),dp)
			end do
		end if
		!
		return
	end function


	
	function kubo_mep_LC(en_k, V_ka, fd_distrib) result(F3)
		!	LOCAL KUBO CONTRIBUTION
		!
		real(dp),			intent(in)		::	en_k(:), fd_distrib(:,:)
		complex(dp),		intent(in)		::	V_ka(:,:,:)
		real(dp),			allocatable		::	F3(:,:,:)
		integer								::	n0, n, n_wf, n_ef, j, k, l, ef_idx, leviC(3,3,3)	
		real(dp) 							::	v_nm_mn_nn(3,3), en_denom
		!
		n_ef	=	size(fd_distrib,1)
		n_wf	=	size(en_k,1)
		allocate(	F3(	3,3,	n_ef))
		!
		call get_levi_civita(leviC)
		!
		F3	= 0.0_dp
		do n0 = 1, n_wf
			!MIXING
			do n = 1, n_wf
			 	if( n/= n0 )	then
			 		!
			 		!	dE	BANDS
			 		en_denom	=	(	en_k(n0) - en_k(n)	)**3		
			 		!
			 		!	TRIPLE PRODUCT
			 		v_nm_mn_nn	=	0.0_dp
			 		do j = 1, 3
			 			do k = 1, 3
			 				do l =1,3
			 					v_nm_mn_nn(:,j)	=	v_nm_mn_nn(:,j) +	 real(	leviC(j,k,l)	,dp) 					&
			 													* real(	V_ka(:,n0,n) * V_ka(k,n,n0) * V_ka(l,n0,n0), dp)
			 				end do
			 			end do
			 		end do
			 		!
			 		!	LOOP FERMI LEVEL
			 		do ef_idx	=	1, n_ef
			 			F3(:,:,ef_idx)	=	F3(:,:,ef_idx)	+ 	v_nm_mn_nn(:,:) * fd_distrib(ef_idx,n0) / en_denom
			 		end do
			 	end if
			 	!
			 	!
			end do
		end do
		!
		return
	end function


	function kubo_mep_IC(en_k, V_ka, fd_distrib) result(F2)
		!	ITINERANT KUBO CONTRIBUTION
		!
		real(dp),			intent(in)		::	en_k(:), fd_distrib(:,:)
		complex(dp),		intent(in)		::	V_ka(:,:,:)
		real(dp),			allocatable		::	F2(:,:,:)
		integer								::	n0, m, n, 				&
												j, k, l, 				&
												n_ef, n_wf, ef_idx,		&
												leviC(3,3,3)
		real(dp) 							::	v_on_nm_mo(3,3), en_denom
		!
		n_ef	=	size(	fd_distrib	,	1)
		n_wf	=	size(	en_k		,	1)
		allocate(	F2(		3,3,	n_ef	))
		!
		F2	= 0.0_dp
		call get_levi_civita(leviC)
		!
		do n0 = 1, n_wf	
			!MIXING
			do m = 1, n_wf
				do n = 1, n_wf
				 	if( n/= n0 .and. m/=n0 )	then
			 			!
			 			!	dE	BANDS
				 		en_denom	=	(en_k(n0) - en_k(m)) 	*	(en_k(n0) - en_k(n))**2	 
				 		!
				 		!	TRIPLE PRODUCT
				 		v_on_nm_mo	=	0.0_dp
				 		do j = 1, 3
				 			do l = 1, 3
				 				do k = 1, 3
				 					v_on_nm_mo(:,j)		=	v_on_nm_mo(:,j)		-	real(leviC(j,k,l),dp)		&
				 									*	real(	V_ka(:,n0,n) * V_ka(k,n,m) * V_ka(l,m,n0)	,dp)
				 				end do
				 			end do
				 		end do
				 		!
				 		!	LOOP FERMI LEVEL
				 		do ef_idx = 1, n_ef
				 			F2(:,:,ef_idx)	=	F2(:,:,ef_idx)	+	v_on_nm_mo(:,:)	* fd_distrib(ef_idx,n0)	/ 	en_denom
				 		end do
				 	end if
				end do
			end do
		end do
		!
		return
	end function




end module kubo_mep