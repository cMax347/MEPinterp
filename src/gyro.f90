module gyro
!!---------------------------------------------------------------------------------------------------------------------
!
!			see arXiv:1710.03204v2
!
!!---------------------------------------------------------------------------------------------------------------------

	use constants,			only:		dp
	use statistics,			only:		fd_stat_deriv
	use matrix_math,		only:		convert_tens_to_vect
	use omp_lib


	implicit none

	public 	::	get_gyro_D, get_gyro_Dw, get_gyro_C, get_gyro_K
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



	function get_gyro_D(V_ka, Om_kab, fd_deriv) result(D_ab)
		complex(dp),					intent(in)			::		V_ka(:,:,:)	
		complex(dp),	allocatable,	intent(in)			::		Om_kab(:,:,:,:)
		real(dp),						intent(in)			::		fd_deriv(:,:)
		complex(dp)											::		om_ca(3), tens(3,3)
		complex(dp),	allocatable							::		D_ab(:,:,:)
		integer												::		n, b, n_ef, n_wf, ef_idx
		!
		n_ef	=	size(fd_deriv,1)
		n_wf	=	size(fd_deriv,2)
		allocate(D_ab(	3,3,	n_ef))
		D_ab	=	0.0_dp
		!
		if(allocated(Om_kab)) then
			!
			!$OMP PARALLEL DO DEFAULT(none) 					& 
			!$OMP PRIVATE(om_ca, tens, b, ef_idx) 				&
			!$OMP SHARED(n_wf, n_ef, V_ka, Om_kab, fd_deriv)	&
			!$OMP REDUCTION(+: D_ab) 
			do n =1 , n_wf
				call	convert_tens_to_vect(Om_kab(:,:,n,n), om_ca(:)	)	
				!
				tens = 0.0_dp
				do b = 1, 3
					tens(:,b)	=	tens(:,b)	+	V_ka(:,n,n) * om_ca(b) 		
				end do
				!
				do ef_idx	=	1 , n_ef
					D_ab(:,:,ef_idx)	=	D_ab(:,:,ef_idx)	+	tens(:,:)	* fd_deriv(ef_idx,n)	
				end do
			end do
			!$OMP END PARALLEL DO
			!
		end if
		!
		return
	end function


	function get_gyro_Dw(en_k, A_ka,hw_lst, fd_deriv)	result(Dw_ab)
		complex(dp),	allocatable,	intent(in)	::	A_ka(:,:,:)
		real(dp),						intent(in)	::	en_k(:), hw_lst(:), fd_deriv(:,:)
		complex(dp),	allocatable					::	Dw_ab(:,:,:,:)
		integer										::	n_hw, n_ef
		!
		!
		n_hw	=	size(hw_lst,1)
		n_ef	=	size(fd_deriv,1)
		!
		allocate(Dw_ab(	3,3,	n_hw,	n_ef))
		Dw_ab	=	0.0_dp
		!
		if(	allocated(A_ka)) then
			if(	size(A_ka,2) /= size(en_k,1))	write(*,*)	"[get_gyro_Dw]: WARNING en_k and A_ka life on different bases"
			!
			!
			!	TODO		!!!!!!!!!!1!
			!
			!
		end if
		!
		return
	end function


	function get_gyro_C(V_ka,  fd_deriv) result(C_ab)
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






end module gyro