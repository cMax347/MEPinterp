module test_matrix_math 
	use parameters,		only:				dp, fp_acc
	use matrix_math,	only:				my_Levi_Civita,                         &
                                            zheevr_wrapper, zheevd_wrapper,         & 
                                            uni_gauge_trafo,                        &
                                            crossP,                                 &
                                            is_equal_vect,                          &
                                            convert_tens_to_vect,                   &
                                            matrix_comm             

    use test_log,		only:				push_to_outFile                         

	implicit none


	private
	public				::				matrix_math_test



	contains

!public:
	logical function matrix_math_test(succ_cnt, nTests)
		integer,			intent(out)			::	succ_cnt, nTests
		logical,			dimension(7)		::	passed
		character(len=80), dimension(7)			::	label
		character(len=10)						:: 	result
		character(len=121)						::	succ_message
		integer									::	test
		!
		nTests	= size(passed)
		succ_cnt= 0
		matrix_math_test	= .true.
		!
		label(		1	)			=	"Levi Civitia operartor"
		passed(		1	)			=	test_levi_civita()
		!-----------------------------------------------------------------
		label(		2	)			=	"equality of 3D vectors"
		passed(		2	)			=	test_vect_equal()
		!-----------------------------------------------------------------
		label(		3	)			=	"tensor to vector conversion"
		passed(		3	)			=	test_tens_vect()
		!-----------------------------------------------------------------
		label(		4	)			=	"Cross product"
		passed(		4	)			=	test_crossp()
		!-----------------------------------------------------------------
		label(		5	)			=	"eigen Solver"
		passed(		5	)			=	test_eig_solver()
		!-----------------------------------------------------------------
		label(		6	)			=	"Gauge rotation"
		passed(		6	)			=	test_gauge_trafo()				
		!-----------------------------------------------------------------			
		label(		7	)			=	"matrix commutator"
		passed(		7	)			=	test_mat_comm()
		!
		do test = 1, nTests
			matrix_math_test	= matrix_math_test .and. passed(test)
			if( passed(test) )	then 
				result	=	"passed"
				succ_cnt=	1 + succ_cnt
			else
				result	=	"not passed"
			end if
			write(succ_message,'(a,a,a,a,a)')	'[matrix_math_test]: ',result,' test of "',trim(label(test)),'"'
			!
			!WRITE TO LOG FILE
			call push_to_outFile(succ_message)
		end do
		!
		write(succ_message,'(a,i1,a,i1,a)')		'[matrix_math_test]: passed ',succ_cnt,' of ',nTests,' tests'
		call push_to_outFile(succ_message)
		call push_to_outFile(" ")
		return
	end function







!private:------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_levi_civita()
		!
		test_levi_civita	= 							(	my_Levi_Civita(1,2,3) ==  1	)	
		test_levi_civita 	= test_levi_civita .and.	(	my_Levi_Civita(2,3,1) ==  1	)
		test_levi_civita	= test_levi_civita .and.	( 	my_Levi_Civita(3,1,2) ==  1	)
		!
		test_levi_civita	= test_levi_civita .and.	( 	my_Levi_Civita(3,2,1) == -1	)
		test_levi_civita	= test_levi_civita .and.	( 	my_Levi_Civita(1,3,2) == -1	)
		test_levi_civita	= test_levi_civita .and.	( 	my_Levi_Civita(2,1,3) == -1	)
		!
		test_levi_civita	= test_levi_civita .and.	( 	my_Levi_Civita(1,1,2) == 0	)
		test_levi_civita	= test_levi_civita .and.	( 	my_Levi_Civita(1,2,2) == 0	)
		test_levi_civita	= test_levi_civita .and.	( 	my_Levi_Civita(2,1,2) == 0	)
		test_levi_civita	= test_levi_civita .and.	( 	my_Levi_Civita(1,1,1) == 0	)
		test_levi_civita	= test_levi_civita .and.	( 	my_Levi_Civita(-1,2,3)== 0	)
		test_levi_civita	= test_levi_civita .and.	( 	my_Levi_Civita(4,1,2) == 0	)
		!
		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_vect_equal()
		real(dp)				::	a(3), b(3), c(3)
		test_vect_equal		=	.false.

		a(:)	= 0.0_dp	
		a(1)	= 1e-14_dp
		!
		b(:)	= 0.0_dp
		b(2)	= -1e-14_dp
		!
		c(:)	= 0.0_dp
		c(1)	= 1e-14_dp
		!
		test_vect_equal	= .not. is_equal_vect(fp_acc, a,b)
		test_vect_equal	=	test_vect_equal .and. is_equal_vect(fp_acc, a,c)


		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_tens_vect()
		test_tens_vect		=	.false.

		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_crossp()
		real(dp)			::	real_a(3), real_b(3), real_c_sol(3), real_c_test(3)
		logical				::	real_test, imag_test
		!
		test_crossp		=	.false.
		real_test		= .false.
		imag_test		= .false.
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
		imag_test	= .true.
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
		test_crossp	= real_test .and. imag_test
		!
		return
	end function	
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_eig_solver()
		complex(dp),	allocatable			::	M_W1(:,:),  M_W2(:,:),	&
												U_l(:,:), U_r(:,:)
		real(dp),		allocatable			::	eigVal(:)
		logical								::	test1, test2
		integer								::	size

		size	= 1000
		!
		allocate(	M_W1(	size, size	)		)
		allocate(	M_W2(	size, size	)		)
		!
		allocate(	U_l(	size, size	)		)
		allocate(	U_r(	size, size	)		)
		allocate(	eigVal(	size		)		)
		!

		! FILL M_W, todo: setup random hermitian matrix
		call get_rand_herm_matrix(M_W1, M_W2, U_l, U_r)									
		!
		call zheevd_wrapper(M_W1, eigVal)
		call zheevd_wrapper(M_W2, eigVal)
		!
		test1	=	equal_mat(M_W1, U_l)
		test2	=	equal_mat(M_W2, U_l)
		!
		test_eig_solver	= test1 .and. test2
		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_gauge_trafo()
		complex(dp),	allocatable			::	M_H(:,:),				& 
												M_W1(:,:),  M_W2(:,:),	&
												M_solv(:,:),			&
												U_l(:,:), U_r(:,:)
		real(dp),		allocatable			::	eigVal(:)
		integer								::	size
		test_gauge_trafo	=	.false.
		!
		size	= 1000
		!
		allocate(	M_H(	size, size	)		)
		allocate(	M_W1(	size, size	)		)
		allocate(	M_W2(	size, size	)		)
		allocate(	M_solv(	size, size	)		)
		!
		allocate(	U_l(	size, size	)		)
		allocate(	U_r(	size, size	)		)
		allocate(	eigVal(	size		)		)
		!
		! FILL M_W, todo: setup random hermitian matrix
		call get_rand_herm_matrix(M_W1, M_W2, U_l, U_r)

		!cpy M_W to U

		!GET U
		M_solv	=	M_W1
		call zheevd_wrapper(M_solv,	eigVal)
		!
		!rotate M_W with U should yield diagonal matrix
		M_H	=	M_W1
		call uni_gauge_trafo(M_solv, M_H)
		!
		!
		!todo: check now if M_H is diagonal with eival on diago


		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_mat_comm()
		!
		!	wikipedia test from 
		!			https://en.wikipedia.org/wiki/Commuting_matrices
		!
		complex(dp),	allocatable		::	I_mat(:,:), M_mat(:,:), cplx_comm(:,:)
		real(dp)						::	A(2,2), B(2,2), real_comm(2,2)
		integer							::	size, m,n, succ_cnt
		real(dp)						::	scale, rand1, rand2
		logical							::	identity_test, wikipedia_comm, wikipedia_non_comm
		!
		!----------------------------------------------------------------------------------------------------------
		!	IDENTITY_TEST:	TEST IF RANDOM MATRIX M COMMUTES WITH IDENTITY
		!
		test_mat_comm	=	.false.
		size			= 	100
		scale			=	1023.34514_dp
		!
		allocate(	I_mat(		size, size	)		)		
		allocate(	M_mat(		size, size	)		)
		allocate(	cplx_comm(	size, size	)		)
		I_mat(:,:)	=	cmplx(0.0_dp, 0.0_dp, dp)
		M_mat(:,:)	=	cmplx(0.0_dp, 0.0_dp, dp)
		!
		!FILL MATRICES
		do n = 1, size
			do m = 1, size
				if(m==n)	I_mat(m,n)	=	cmplx(1.0,0.0_dp, dp)
				call random_number(	rand1	)			
				call random_number(	rand2	)
				M_mat(m,n)	=	cmplx( 		scale*(	rand1-.5_dp)	,	scale*(	rand2-.5_dp),	dp)
			end do
		end do
		!
		!GET COMMUTATOR
		call matrix_comm(I_mat, M_mat, cplx_comm)
		!
		!TEST IF COMMUTATOR IS ZERO
		succ_cnt = 0
		do n = 1, size
			do m = 1, size
				if(	 abs(	cplx_comm(m,n)	) < fp_acc		)	succ_cnt = succ_cnt + 1
			end do
		end do
		identity_test	=	( succ_cnt == size**2 )	
		!----------------------------------------------------------------------------------------------------------
		!	wikipedia_comm:	COMMUTING MATRICES TEST
		!
		A(:,:)	=	0.0_dp
		A(1,2)	=	1.0_dp
		B(:,:)	=	0.0_dp
		B(1,1)	=	1.0_dp
		B(2,2)	=	1.0_dp
		!
		call matrix_comm(A,B,	real_comm )
		wikipedia_comm	=							(	abs(real_comm(1,1))	< fp_acc	)
		wikipedia_comm	=	wikipedia_comm	.and.	(	abs(real_comm(1,2))	< fp_acc	)
		wikipedia_comm	=	wikipedia_comm	.and.	(	abs(real_comm(2,1))	< fp_acc	)
		wikipedia_comm	=	wikipedia_comm	.and.	(	abs(real_comm(2,2))	< fp_acc	)
		!----------------------------------------------------------------------------------------------------------
		!	wikipedia_comm:	NON COMMUTING MATRICES TEST
		!
		A(1,1)	=	1.0_dp
		A(1,2)	=	2.0_dp
		A(2,1)	=	0.0_dp
		A(2,2)	=	3.0_dp

		B(:,:)	=	1.0_dp
		B(2,1)	=	0.0_dp
		!
		call matrix_comm(A,B,	real_comm )
		wikipedia_non_comm	=								(	abs(real_comm(1,1))	< fp_acc	)
		wikipedia_non_comm	=	wikipedia_non_comm	.and.	(	abs(real_comm(1,2))	< fp_acc	)
		wikipedia_non_comm	=	wikipedia_non_comm	.and.	(	abs(real_comm(2,1))	< fp_acc	)
		wikipedia_non_comm	=	wikipedia_non_comm	.and.	(	abs(real_comm(2,2))	< fp_acc	)
		!
		wikipedia_non_comm	=	.not. wikipedia_non_comm
		!----------------------------------------------------------------------------------------------------------
		!
		test_mat_comm	=	identity_test .and. wikipedia_comm .and. wikipedia_non_comm
		!
		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------





!helpers:

	subroutine get_rand_herm_matrix(Mat_a, Mat_b, left_eig, right_eig)
		!	
		!	https://software.intel.com/en-us/mkl-developer-reference-fortran-latm6#045E00E9-CEEC-40C3-8283-0CE8166BF2B9
		!
		complex(dp), intent(out)	::	Mat_a(:,:), Mat_b(:,:), left_eig(:,:), right_eig(:,:)
		integer						::	typ, n, lda, ldx, ldy
		complex(dp)					::	alpha, beta, wx, wy		
		complex(dp),	allocatable	::	cond_eig_val(:), cond_eig_vec(:)
		!
		typ		= 2
		n		= size(Mat_a		, 1)
		lda		= size(Mat_a		, 1)
		ldx		= size(right_eig	, 1)
		ldy		= size(left_eig		, 1)
		!
		allocate(		cond_eig_val( n )		)		
		allocate(		cond_eig_vec( n )		)
		!
		alpha	=	cmplx(1.3_dp		, 	-0.4_dp		,dp)
		beta	=	cmplx(-3.1_dp		, 	2.5_dp		,dp)
		wx		=	cmplx(1.3_dp		, 	-0.0013_dp	,dp)
		wy		=	cmplx(-0.0001_dp	, 	-2.567_dp	,dp)
		!
		!call zlatm6( type, n, a, lda, b, x, ldx, y, ldy, alpha, beta, wx, wy, s, dif )
		call zlatm6( typ, n, Mat_a, lda, Mat_b, right_eig, ldx, left_eig, ldy, alpha, beta, wx, wy, cond_eig_val, cond_eig_vec)
		!
		return
	end subroutine	



	logical function	equal_mat(A, B)
		complex(dp),	intent(in)		::	A(:,:), B(:,:)
		real(dp)						::	delta
		integer							::	n, m
		!
		equal_mat	=	(	size(A,1) /= size(B,1) .or. size(A,2) /= size(B,2) )
		!
		if( equal_mat	) then

			loop_out: do n = 1, size(A,2)
				loop_in: do m = 1, size(A,1)
					delta	=	abs(	A(m,n) - B(m,n)	)
					!
					equal_mat	= equal_mat .and. (delta < fp_acc)
					!
					if(.not. equal_mat) exit loop_out
				end do loop_in
			end do loop_out
		end if
		!
		return
	end function





































end module               