module statistics
	use constants, 		only:			dp, kBoltz_Eh_K

	implicit none

	public				::	fd_stat, fd_get_N_el
	private

contains


	real(dp) pure function fd_stat(e_band, e_fermi,	T_kelvin) 
		real(dp), 		intent(in)			::	e_band, e_fermi, T_kelvin
		real(dp)							::	T_smear
		!
		T_smear			=	kBoltz_Eh_K		*	T_kelvin
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

end module statistics