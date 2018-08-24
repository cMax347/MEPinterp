program test_crossp
	use helpers,			only:			my_exit,				&
											init_outFile, 			&
											push_to_outFile
	use parameters,			only:			dp, fp_acc,				& 
											crossP


	implicit none


	logical			::	passed, real_test, imag_test
	integer			::	cnt
	!
	!
	cnt				=	0
	!	TESTS
	real_test		=	d_test_crossp()
	imag_test		= 	z_test_crossp()
	!
	!
	!	LOG
	if( real_test	)	then	
							call push_to_outFile('[test_crossp]: PASSED real cross product')
	else
							call push_to_outFile("[test_crossp]: FAILED real cross product")
	end if
	if( imag_test	)	then
							call push_to_outFile('[test_crossp]: PASSED imag cross product')
	else
							call push_to_outFile("[test_crossp]: FAILED	complex cross product")
	end if
	call push_to_outFile("------------------------------------------------------")
	call push_to_outFile("")
	!
	!
	!	EXIT
	passed = real_test .and. imag_test
	call my_exit(passed)


contains


	logical function d_test_crossp()
		real(dp)		::	real_a(3), real_b(3), real_c_sol(3), real_c_test(3)
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
		d_test_crossp	= 	norm2(real_c_test - real_c_sol) < fp_acc 
		!
		return
	end function

	logical function z_test_crossp()
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
		!z_test_crossp	=	abs(	sum(cplx_delta)		) < fp_prec
		!
		z_test_crossp	=	.true.
		!
		return
	end function


	
end program


