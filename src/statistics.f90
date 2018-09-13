module statistics
	use constants, 		only:			dp, kBoltz_Eh_K



	implicit none


	public				::	fermi_dirac

	private



contains


	real(dp) pure function fermi_dirac(e_band, e_fermi,	T_kelvin) 
		real(dp), 		intent(in)			::	e_band, e_fermi, T_kelvin
		real(dp)							::	e_temp
		!
		e_temp	=	kBoltz_Eh_K	*	T_kelvin
		!
		fermi_dirac = 1.0_dp	/	(		1.0_dp	+	exp(	(e_band	- e_fermi)	/	(e_temp))				)
		!
	end function



end module statistics