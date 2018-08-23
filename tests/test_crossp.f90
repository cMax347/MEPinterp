program test_crossp
	use test_log,			only:			init_outFile, 			&
											push_to_outFile
	use parameters,			only:			dp, fp_acc,				& 
											crossP


	implicit none


	real(dp)		::	real_a(3), real_b(3), real_c_sol(3), real_c_test(3)
	logical			::	passed, real_test, imag_test
	passed			=	.false.
	real_test		= 	.false.
	imag_test		= 	.false.
	!
	!REAL(DP) TEST
	real_a(1)		= 	-1.4_dp
	real_a(2)		=	3.5_dp
	real_a(3)		=	2.6_dp
	!
	real_b(1)		=	0.9_dp	
	real_b(2)		=	-14.5_dp	
	real_b(3)		=	0.001_dp		
	!
	real_c_sol(1)	=	37.703500_dp
	real_c_sol(2)	=	2.341400_dp
	real_c_sol(3)	=	17.1500_dp
	!
	real_c_test		=	crossP(real_a,real_b)
	real_test		= 	norm2(real_c_test - real_c_sol) < fp_acc 
	!
	!COMPLEX(DP) TEST
	imag_test		= .true.
	!cplx_a(1)		=	cmplx(-1.4_dp, 1.0_dp,dp)
	!cplx_a(2)		=	cmplx(3.5_dp,-0.5_dp,dp)
	!cplx_a(3)		=	cmplx(2.6_dp,0.0_dp,dp)
	!!
	!cplx_b(1)		=	cmplx(0.9_dp,0.0_dp,dp)
	!cplx_b(2)		=	cmplx(-14.5_dp,0.0_dp,dp)
	!cplx_b(3)		=	cmplx(0.001_dp,-0.001_dp,dp)
	!!
	!cplx_c_sol(1)	=	cmplx(37.7030_dp,0.004_dp,dp)
	!cplx_c_sol(2)	=	cmplx( 2.3404_dp,0.0024_dp,dp)
	!cplx_c_sol(3)	=	cmplx(17.15_dp,14.05_dp,dp)
	!!
	!cplx_c_test		=	crossP(cplx_a, cplx_b)
	!cplx_delta		=	cplx_c_test - cplx_c_sol
	!imag_test		=	abs(	sum(cplx_delta)		) < fp_prec
	!
	if(.not. real_test	)	call push_to_outFile("[test_crossp]: FAILED real    cross product")
	if(.not. imag_test	)	call push_to_outFile("[test_crossp]: FAILED	complex cross product")
	!
	passed = real_test .and. imag_test



	if(passed)	then
		call exit(0)
	else
		call exit(1)
	end if




	
end program


