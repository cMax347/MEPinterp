module constants
	implicit none

	
	!FLOATING POINT
	integer, 		parameter 	:: 	dp 				= kind(0.d0)
	real(dp),		parameter	:: 	fp_acc			= 1e-14_dp

	!MPI
	integer						::	mpi_id, mpi_root_id, mpi_nProcs, ierr

	!MATHEMATICAL
	real(dp), 		parameter 	::	PI_dp 			= 4 * atan (1.0_dp)
	complex(dp),	parameter 	::	i_dp 			= cmplx(0.0_dp, 1.0_dp, dp)


	!PHYSICAL
	real(dp),		parameter 	::	aUtoAngstrm 	= 0.52917721092_dp, &
									aUtoEv	 		= 27.211385_dp
	real(dp),		parameter	::	aUtoTesla		= 235051.76_dp
	real(dp),		parameter	::	speedOfLight	= 137.035999_dp !in atomic units
	



end module constants