module gyro
!!---------------------------------------------------------------------------------------------------------------------
!
!			see arXiv:1710.03204v2
!
!!---------------------------------------------------------------------------------------------------------------------

	use constants,			only:		dp
	use statistics,			only:		fd_stat_deriv
	use matrix_math,		only:		convert_tens_to_vect



	implicit none

	public 	::	get_gyro_D, get_gyro_Dw, get_gyro_C, get_gyro_K
	private


	
contains

	pure function get_gyro_K(en, v_a, m_a, e_fermi, T_kelvin) result(K_ab)
		real(dp),		intent(in)			::		en(:), m_a(:,:),	e_fermi, T_kelvin
		complex(dp),	intent(in)			::		v_a(:,:,:)
		complex(dp)							::		K_ab(3,3)
		integer								::		n, b
		!
		K_ab	=	0.0_dp
		!
		do n = 1, size(v_a,3)
			!
			!
			do b = 1, 3
				K_ab(:,b)	=	K_ab(:,b)	+	v_a(:,n,n) * m_a(b,n) 		* fd_stat_deriv(en(n), e_fermi, T_kelvin) 	
			end do
			!
			!
		end do
		!
		return
	end function



	pure function get_gyro_D(en, v_a, om_cab, e_fermi, T_kelvin) result(D_ab)
		real(dp),						intent(in)			::		en(:), e_fermi, T_kelvin
		complex(dp),					intent(in)			::		v_a(:,:,:)	
		complex(dp),	allocatable,	intent(in)			::		Om_cab(:,:,:,:)
		complex(dp)											::		om_ca(3)
		complex(dp)											::		D_ab(3,3)
		integer												::		n, b
		!
		D_ab	=	0.0_dp
		!
		if(allocated(Om_kab)) then
			do n =1 , size(om_cab,	4)
				call	convert_tens_to_vect(om_cab(:,:,n,n), om_ca(:)	)	
				!
				!
				do b = 1, 3
					D_ab(:,b)	=	D_ab(:,b)	+	v_a(:,n,n) * om_ca(b) 		* fd_stat_deriv(en(n), e_fermi, T_kelvin)		
				end do
				!
				!
			end do
		end if
		!
		return
	end function


	pure function get_gyro_Dw()	result(Dw_ab)
		complex(dp)							::		Dw_ab(3,3)
		!
		Dw_ab	=	0.0_dp
		!
		return
	end function


	pure function get_gyro_C(en, v_a,  e_fermi, T_kelvin) result(C_ab)
		real(dp),		intent(in)			::		en(:), e_fermi, T_kelvin
		complex(dp),	intent(in)			::		v_a(:,:,:)
		complex(dp)							::		C_ab(3,3)
		integer								::		n, b
		!
		C_ab	=	0.0_dp
		!
		do n = 1, size(v_a, 3)
			!
			!
			do b = 1, 3
				C_ab(:,b)	=	C_ab(:,b)	+	v_a(:,n,n) * v_a(b,n,n)			* fd_stat_deriv(en(n), e_fermi, T_kelvin)
			end do
			!
			!
		end do
		!
		return
	end function






end module gyro