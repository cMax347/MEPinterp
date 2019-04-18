module constants
	implicit none
	save
	!
	!FLOATING POINT
	integer, 		parameter 	:: 	dp 				= kind(0.d0)
	real(dp),		parameter	:: 	fp_acc			= 1e-14_dp
	real(dp),		parameter	::	real_0			= 0.0_dp
	complex(dp),	parameter	:: 	cmplx_0			= cmplx(real_0,real_0,dp) 
	!		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	!MATHEMATICAL
	real(dp), 		parameter 	::	PI_dp 			= 4 * atan (1.0_dp)
	complex(dp),	parameter 	::	i_dp 			= cmplx(0.0_dp, 1.0_dp, dp)
	!		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	!PHYSICAL
	real(dp),		parameter 	::	aUtoAngstrm 	= 0.52917721092_dp, &			!	Bohr 			->	Angstroem
									aUtoEv	 		= 27.211385_dp					!	Hartree			->	eV
	real(dp),		parameter	::	aUtoTesla		= 235051.76_dp					!	mag.Field(au,si)->	Tesla
	real(dp),		parameter	::	speedOfLight	= 137.035999_dp 				!	in atomic units
	real(dp),		parameter	::	kBoltz_Eh_K		= 8.6173303_dp * 1e-5_dp / aUtoEv 		!kBoltz from wiki: kB	=	8.6173303(50)×10−5	eV⋅K−1
	real(dp),		parameter	::	T_room			= 300.0_dp
	!		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end module constants