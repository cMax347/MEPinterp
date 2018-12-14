module statistics
	use constants, 		only:			dp, kBoltz_Eh_K

	implicit none

	public					::			fd_stat,		fd_stat_deriv, 			&
										fd_get_N_el, fd_count_el
	private

contains


	pure subroutine fd_count_el(en_k, eFermi, T_kelvin, el_count, sum_loc, n_el_min, n_el_max)
		real(dp),		intent(in)					::		en_k(:), eFermi, T_kelvin
		real(dp),		intent(out)					::		el_count
		real(dp),		intent(inout)				::		sum_loc, n_el_min, n_el_max
		!
		el_count			=	fd_get_N_el(en_k, eFermi, T_kelvin)
		sum_loc				=	sum_loc		+	el_count
		!		!
		if(	el_count	<	n_el_min	)		n_el_min	= 	el_count
		if( el_count	>	n_el_max	)		n_el_max	=	el_count
		!
		return
	end subroutine



	real(dp) pure function fd_stat(e_band, e_fermi,	T_kelvin) 
		real(dp), 		intent(in)			::	e_band, e_fermi, T_kelvin
		real(dp)							::	T_smear
		!
		T_smear			=	kBoltz_Eh_K		*	T_kelvin
		!
		if(T_smear < 1e-6_dp ) T_smear = 1e-6_dp

		!
		fd_stat		 	= 	1.0_dp	/	(		1.0_dp	+	exp(	(e_band	- e_fermi)	/	(T_smear))				)
		!
		return
	end function




	real(dp) pure function fd_get_N_el(en_k, e_fermi, T_kelvin)
		real(dp),		intent(in)			::	en_k(:), e_fermi, T_kelvin
		integer								::	n
		!
		fd_get_N_el	=	0.0_dp
		!
		do n = 1, size(en_k)
			fd_get_N_el	=	fd_get_N_el	+	fd_stat(en_k(n),e_fermi, T_kelvin)
		end do
		!
		return
	end function



	real(dp) pure function fd_stat_deriv(e_band, e_fermi, T_kelvin)
		real(dp), 		intent(in)			::	e_band, e_fermi, T_kelvin
		real(dp)							::	T_smear, x
		!
		!	w90git:
		!if (abs (x) .le.36.0) then
        ! 	utility_w0gauss = 1.00_dp / (2.00_dp + exp ( - x) + exp ( + x) )
        !  	! in order to avoid problems for large values of x in the e
       	!else
        !	utility_w0gauss = 0.0_dp
       	!endif
       	!
       	!
       	T_smear			=	kBoltz_Eh_K		*	T_kelvin
       	x				=	(e_band	- e_fermi) / T_smear
       	!
       	!
       	fd_stat_deriv	=	1.00_dp 	/ 		(	2.00_dp + exp( x ) + exp( -x ) 			)

		return
	end function	

end module statistics