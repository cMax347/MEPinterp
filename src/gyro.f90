module gyro
!!---------------------------------------------------------------------------------------------------------------------
!
!			see arXiv:1710.03204v2
!	or
!			Souza et al., PRB 97 035158 (2018)
!!---------------------------------------------------------------------------------------------------------------------

	use constants,			only:		dp
	use statistics,			only:		fd_stat_deriv
	use matrix_math,		only:		convert_tens_to_vect, &
										crossP
	use wann_interp,		only:		get_kubo_curv
	use omp_lib


	implicit none

	public 						::		get_gyro_D, get_gyro_Dw,	& 
										get_gyro_C, 				&
										get_gyro_K
	private


	
contains

	function get_gyro_K(V_ka, m_a, fd_distrib) result(K_ab)
		real(dp),		intent(in)			::		m_a(:,:),	fd_distrib(:,:)
		complex(dp),	intent(in)			::		V_ka(:,:,:)
		complex(dp),	allocatable			::		K_ab(:,:,:)
		complex(dp)							::		tens(3,3)
		integer								::		n, b, n_wf, n_ef, ef_idx
		!
		!	m_a is magnetization 
		!
		!	can we use equation (15) for the orbital magnetization
		!
		n_wf	=	size(fd_distrib,2)
		n_ef	=	size(fd_distrib,1)
		allocate(	K_ab(	3,3,	n_ef	))
		K_ab	=	0.0_dp
		!
		!$OMP PARALLEL DO  DEFAULT(none) 				&
		!$OMP PRIVATE(tens, b, ef_idx)					& 
		!$OMP SHARED(n_wf, n_ef, V_ka, m_a, fd_distrib)	&
		!$OMP REDUCTION(+: K_ab)
		do n = 1, n_wf
			!
			!
			tens = 0.0_dp
			do b = 1, 3
				tens(:,b)	=	tens(:,b)	+	V_ka(:,n,n) * m_a(b,n) 		
			end do
			!
			do ef_idx	= 1, n_ef
				K_ab(:,:,ef_idx)	=	K_ab(:,:,ef_idx)	+	tens(:,:)	* fd_distrib(ef_idx,n)
			end do
			!
			!
		end do
		!$OMP END PARALLEL DO
		!
		return
	end function



	function get_gyro_D(en_k, V_ka, Om_kab, fd_deriv) result(D_ab)
		!
		!	calculates Berry curvature dipole matrix D	
		!
		!		->	D is traceless
		!		->	D is dimensionless
		!
		!	see eq.(2) of PRB 97 035158 (2018)
		!
		real(dp),						intent(in)			::		en_k(:), fd_deriv(:,:)
		complex(dp),					intent(in)			::		V_ka(:,:,:)	
		complex(dp),	allocatable,	intent(in)			::		Om_kab(:,:,:,:)
		complex(dp),	allocatable							::		kubo_curv(:,:,:)		
		complex(dp)											::		om_ca(3)
		complex(dp)											::		tens(3,3)
		complex(dp),	allocatable							::		D_ab(:,:,:)
		integer												::		n, b, n_ef, n_wf, ef_idx
		!
		n_ef	=	size(fd_deriv,1)
		n_wf	=	size(fd_deriv,2)
		allocate(	D_ab(	3,3, n_ef))
		D_ab	=	0.0_dp
		!
		if(	.not.	allocated(Om_kab) )	then
			!	USE KUBO TO ESTIMATE CURV
			call get_kubo_curv(	en_k, V_ka, kubo_curv	)
			!
		else
			!	USE GIVEN CURVATURE
			allocate(kubo_curv(3,3,n_wf))
			do n = 1, n_wf
				kubo_curv(:,:,n)	=	Om_kab(:,:,n,n)
			end do
		end if
		!
		!$OMP PARALLEL DO DEFAULT(none) 					& 
		!$OMP PRIVATE(om_ca, tens, b, ef_idx) 				&
		!$OMP SHARED(n_wf, n_ef, V_ka, kubo_curv, fd_deriv)	&
		!$OMP REDUCTION(+: D_ab) 
		do n =1 , n_wf
			tens	=	kubo_curv(:,:,n)
			call	convert_tens_to_vect(tens(:,:), om_ca(:)	)	
			!
			tens = 0.0_dp
			do b = 1, 3
				tens(:,b)	=	tens(:,b)	+	V_ka(:,n,n) * om_ca(b) 		
			end do
			!
			do ef_idx	=	1 , n_ef
				D_ab(:,:,ef_idx)	=	D_ab(:,:,ef_idx)	-	tens(:,:)	* fd_deriv(ef_idx,n)
			end do
		end do
		!$OMP END PARALLEL DO
		!!
		!else
		!	!
		!	!$OMP PARALLEL DO DEFAULT(none) 					& 
		!	!$OMP PRIVATE(om_ca, tens, b, ef_idx) 				&
		!	!$OMP SHARED(n_wf, n_ef, V_ka, Om_kab, fd_deriv)	&
		!	!$OMP REDUCTION(+: D_ab) 
		!	do n =1 , n_wf
		!		call	convert_tens_to_vect(Om_kab(:,:,n,n), om_ca(:)	)	
		!		!
		!		tens = 0.0_dp
		!		do b = 1, 3
		!			tens(:,b)	=	tens(:,b)	+	V_ka(:,n,n) * om_ca(b) 		
		!		end do
		!		!
		!		do ef_idx	=	1 , n_ef
		!			D_ab(:,:,ef_idx)	=	D_ab(:,:,ef_idx)	-	tens(:,:)	* fd_deriv(ef_idx,n)
		!		end do
		!	end do
		!	!$OMP END PARALLEL DO
		!	!
		!end if
		!
		return
	end function


	function get_gyro_Dw(en_k, V_ka, A_ka,hw_lst, fd_deriv)	result(Dw_ab)
		!
		!	calculates the finite-frequency generalization of Berry curvature dipole
		!
		!	see eq.(12) of PRB 97 035158 (2018)
		!
		complex(dp),	allocatable,	intent(in)	::	V_ka(:,:,:), A_ka(:,:,:)
		real(dp),						intent(in)	::	en_k(:), hw_lst(:), fd_deriv(:,:)
		complex(dp),	allocatable					::	Dw_ab(:,:,:,:), tens(:,:,:)
		real(dp),		allocatable					::	Om_hw(:,:,:)
		integer										::	n_wf, n_hw, n_ef, n, hw, b, ef_idx
		!
		!
		n_wf	=	size(en_k,		1)
		n_hw	=	size(hw_lst,	1)
		n_ef	=	size(fd_deriv,	1)
		!
		allocate(	Om_hw(		  3,	n_hw,			n_wf	))
		allocate(	tens(		3,3,	n_hw					))
		allocate(	Dw_ab(		3,3,	n_hw,	n_ef			))
		!
		Dw_ab	=	0.0_dp
		!
		if(	allocated(A_ka)) then
			!	PRECALC CURVATURE( FREQUENCY W)
			call get_curvTilde_hw(	en_k, A_ka, hw_lst, Om_hw)
			!
			!	LOOP BANDS
			do n =1 , n_wf
				!	ITERATE REAL SPACE DIRECTIONS				
				tens 	=	0.0_dp
				do hw = 1, n_hw
					do b = 1, 3
						tens(:,b,hw)	=	tens(:,b,hw)	+	V_ka(:,n,n) * Om_hw(b,hw,n)
					end do
				end do
				!
				!	ITERATE FERMI LEVEL
				do ef_idx = 1, n_ef
					Dw_ab(:,:,:,ef_idx)	=	Dw_ab(:,:,:,ef_idx)	-	tens(:,:,:)	*	 fd_deriv(ef_idx, n)
				end do
			end do 
			!
		end if
		!
		return
	end function


	function get_gyro_C(V_ka,  fd_deriv) result(C_ab)
		!
		!	calculates the Ohmic conductivity in constant relaxation time approximation
		!	units of surface current density (A/cm)
		!
		real(dp),		intent(in)			::		fd_deriv(:,:)
		complex(dp),	intent(in)			::		V_ka(:,:,:)
		complex(dp),	allocatable			::		C_ab(:,:,:)
		complex(dp)							::		tens(3,3)

		integer								::		n, b, n_ef, n_wf, ef_idx
		!
		n_ef	=	size(fd_deriv,1)
		n_wf	=	size(fd_deriv,2)
		allocate(	C_ab(	3,3,	n_ef	))
		C_ab	=	0.0_dp
		!
		!$OMP PARALLEL DO DEFAULT(none)		 	&
		!$OMP PRIVATE(b, tens, ef_idx) 			&
		!$OMP SHARED(n_wf, n_ef, V_ka, fd_deriv)	&
		!$OMP REDUCTION(+: C_ab)
		do n = 1, n_wf
			!
			!
			tens = 0.0_dp
			do b = 1, 3
				tens(:,b)	=	tens(:,b)	+	V_ka(:,n,n) * V_ka(b,n,n)			
			end do
			!
			do ef_idx = 1, n_ef
				C_ab(:,:,ef_idx)	=	C_ab(:,:,ef_idx)	+	tens(:,:)	* fd_deriv(ef_idx, n)
			end do
			!
		end do
		!$OMP END PARALLEL DO
		!
		return
	end function


!	private:

	subroutine get_curvTilde_hw(en_k, A_ka, hw_lst, Om_tilde)
		!
		!	see eq.(C20) of PRB 97 035158 (2018)
		!
		real(dp),		intent(in)				::	en_k(:), hw_lst(:)
		complex(dp),	intent(in)				::	A_ka(:,:,:)
		real(dp)								::	om_tilde(:,:,:)	
		real(dp)								::	ww_kmn
		real(dp)								::	om_vect(3)
		integer									::	n, m, n_wf, hw, n_hw
		!
		n_wf	=	size(en_k,	1)
		n_hw	=	size(hw_lst,1)
		!
		do n = 1, n_wf
			do m = 1, n_wf
				if(	m==n )	cycle
				!
				ww_kmn		=	(	en_k( m) -	en_k( n)	)**2
				om_vect(:)	=	aimag(	crossP(	A_ka(:,n,m),	A_ka(:,m,n)	)	)	
				!
				!
				om_vect = 	om_vect * ww_kmn
				do hw = 1, n_hw
					om_tilde(:,hw,n)	=	om_tilde(:, hw, n)	- 	om_vect(:)	 /	(	ww_kmn - hw_lst(hw)	)
				end do
			end do
		end do
		!
		return
	end subroutine


end module gyro
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~
!~