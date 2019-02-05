program test_matrix_math 
	use constants,		only:				dp, fp_acc
	use matrix_math,	only:				get_levi_civita, 	&
											crossP,				&
											is_equal_vect,		&
											is_equal_mat,		&
											zheevd_wrapper,		&
											blas_matmul,		&
											matrix_comm

    use helpers,		only:				my_exit,									&
    										push_to_outFile, write_test_results, 		&
    										random_matrix, random_vector                     

	implicit none


	logical		::	all_passed
	integer		::	smpl_size, succ_cnt, nTests
	real(dp)	::	acc_limit=1e-11_dp
	logical,		   allocatable		::	passed(:)
	character(len=80), allocatable		::	label(:)
	!
	all_passed	= 	.false.
	nTests		= 	7
	succ_cnt	= 	0
	smpl_size	=	500	
	allocate(	passed(nTests)	)		
	allocate(	label(nTests)	)		
	!
	!	PERFORM TESTS
	!
	label(		1	)			=	"Levi Civitia operartor"
	passed(		1	)			=	test_levi_civita()
	!-----------------------------------------------------------------
	label(		2	)			=	"3D cross product"
	passed(		2	)			=	test_cross_product()
	!-----------------------------------------------------------------
	label(		3	)			=	"equality of  vectors"
	passed(		3	)			=	test_vect_equal()
	!-----------------------------------------------------------------
	label(		4	)			=	"equality of  matrices"
	passed(		4	)			=	test_mat_equal()
	!-----------------------------------------------------------------
	label(		5	)			=	"Gauge rotation"
	passed(		5	)			=	test_gauge_trafo()				
	!-----------------------------------------------------------------			
	label(		6	)			=	"matmul with Blas"
	passed(		6	)			=	test_blas_matmul()
	!-----------------------------------------------------------------			
	label(		7	)			=	"matrix commutator"
	passed(		7	)			=	test_mat_comm()
	!
	!
	!	WRITE RESULTS TO LOG
	!
	call write_test_results(	'matrix_math_test', passed, label, succ_cnt, nTests	)
	call push_to_outFile("------------------------------------------------------")
	call push_to_outFile("")
	!
	!	EXIT
	!
	all_passed	=	( succ_cnt	==	nTests)
	call my_exit(all_passed)


contains







!private:------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_levi_civita()
		integer				::	my_Levi_Civita(3,3,3)
		!
		call get_levi_civita(my_Levi_Civita)
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
		!
		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!

	logical function test_cross_product()
		logical				::	re_crossp, im_crossp
		!
		re_crossp			=	d_test_crossp()
		im_crossp			=	z_test_crossp()
		!
		test_cross_product	=	re_crossp .and. im_crossp
		return
	end function

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


!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_vect_equal()
		real(dp), 	allocatable				::	a(:), b(:), a_cpy(:)
		!integer								::	smpl_size
		!
		!smpl_size	= 1000
		allocate( a(		smpl_size)	)
		allocate( b(		smpl_size)	)		
		allocate( a_cpy(	smpl_size)	)
		!
		call random_vector(a)
		call random_vector(b)
		b(1)		=	a(1)-acc_limit
		a_cpy(:)	= 	a(:)

		!
		test_vect_equal	=  .not. is_equal_vect(fp_acc, a,b)
		test_vect_equal	=	test_vect_equal .and. is_equal_vect(fp_acc, a,a_cpy)
		!
		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!	
	logical function test_mat_equal()
		logical							::	real_test, imag_test
		!
		real_test	=	d_test_mat_equal()
		imag_test	=	z_test_mat_equal()
		!
		if(.not.	real_test)	call push_to_outFile('[test_mat_equal]:  FAILED real test')		
		if(.not.	real_test)	call push_to_outFile('[test_mat_equal]:  FAILED cplx test')
		!
		test_mat_equal	=	real_test .and. imag_test
		!
		return
	end function


	!-------
	logical function d_test_mat_equal()
		real(dp),		allocatable		::	A(:,:), B(:,:), A_cpy(:,:)
		real(dp)						::	acc
		character(len=120)				:: 	msg
		!
		allocate(		A(		smpl_size,	smpl_size)		)
		allocate(		A_cpy(	smpl_size,	smpl_size)		)
		allocate(		B(	smpl_size,	smpl_size)		)
		!
		d_test_mat_equal	= .false.
		!
		call random_matrix(A)
		A_cpy	=	A
		!	now get B matrix
		call random_matrix(B)
		!		(make explicitly sure that it is not equal A)
		B(smpl_size,smpl_size)	=	A(smpl_size,smpl_size)	+ 1.0_dp
		!
		!
		acc	= fp_acc
		do while(.not. d_test_mat_equal .and. acc< acc_limit)
			d_test_mat_equal	= 	is_equal_mat(acc,A, A_cpy) .and. (	.not. is_equal_mat(acc,A,B)	)
			acc					=	acc * 10.0_dp
		end do
		!
		if( 		d_test_mat_equal)	write(msg,*) '[d_test_mat_equal]: test was successful with accuracy ',acc
		if( .not. 	d_test_mat_equal)	write(msg,*) '[d_test_mat_equal]: FAILED with final accuracy ',acc
		call push_to_outFile(msg)
		!
		return
	end function

	logical function z_test_mat_equal()
		complex(dp),		allocatable		::	A(:,:), B(:,:), A_cpy(:,:)
		real(dp)							::	acc
		character(len=120)					:: 	msg
		!
		z_test_mat_equal	= .false.
		!
		allocate(		A(		smpl_size,	smpl_size)		)
		allocate(		A_cpy(	smpl_size,	smpl_size)		)
		allocate(		B(	smpl_size,	smpl_size)		)
		!
		call random_matrix(A)
		A_cpy	=	A
		!	now get B matrix
		call random_matrix(B)
		!		(make explicitly sure that it is not equal A)
		B(smpl_size,smpl_size)	=	A(smpl_size,smpl_size)	+ 1.0_dp
		!
		!
		acc	= fp_acc
		do while(.not. z_test_mat_equal .and. acc< acc_limit)
			z_test_mat_equal	= is_equal_mat(acc,A, A_cpy) .and. .not. is_equal_mat(acc,A,B)
			acc					=	acc * 10.0_dp
		end do
		!
		if( 		z_test_mat_equal)	write(msg,*) '[z_test_mat_equal]: test was successful with accuracy ',acc
		if( .not. 	z_test_mat_equal)	write(msg,*) '[z_test_mat_equal]: FAILED with final accuracy ',acc
		call push_to_outFile(msg)
		!
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
	logical function test_gauge_trafo()
		complex(dp),	allocatable			::	M_W(:,:),				&
												U(:,:), tmp(:,:)
		real(dp),		allocatable			::	eigVal(:)
		real(dp)							::	acc
		character(len=120)					::	msg
		!
		test_gauge_trafo	=	.false.
		!
		allocate(	M_W(	smpl_size, smpl_size	)		)
		allocate(	tmp(	smpl_size, smpl_size	)		)
		allocate(	U(		smpl_size, smpl_size	)		)
		allocate(	eigVal(	smpl_size		)		)
		!
		call random_matrix(M_W)
		!GET EIGVAL & EIGVEC
		U	=	M_W
		call zheevd_wrapper(U,	eigVal)

		!rotate M_W with U should yield diagonal matrix (matrix in H gauge)
		call blas_matmul(	conjg(transpose(U)),	M_W,		tmp)
		call blas_matmul(	tmp,					U,			M_W)
		!
		!
		!todo: check now if M_H is diagonal with eival on diago
		acc = fp_acc
		do while( .not. test_gauge_trafo .and. (acc	<	acc_limit)	)
			test_gauge_trafo	=	test_gauge_with_acc(acc, M_W, eigVal)
			acc 				=	acc * 10.0_dp
		end do
		!
		if( test_gauge_trafo)	 then
			write(msg,*)	'[test_gauge_trafo]:	PASSED gauge rotation with accuracy', acc
			call push_to_outFile(msg)
		else
			write(msg,*)	'[test_gauge_trafo]:	FAILED gauge rotation with final accuracy', acc
			call push_to_outFile(msg)
		end if


		return
	end function

	logical function test_gauge_with_acc(acc, M_W, eigVal)
		real(dp),		intent(in)			::	acc, eigVal(:)
		complex(dp),	intent(in)			::	M_W(:,:)
		integer								::	n,m
		!
		test_gauge_with_acc	= .true.
		outer_loop:	do m = 1, size(M_W,2)
			do n = 1, size(M_W,1)
				if(	n==m	) then
					test_gauge_with_acc =	test_gauge_with_acc .and. 		(	abs(	M_W(n,n)-eigVal(n) 	) < acc	)		
				else
					test_gauge_with_acc =	test_gauge_with_acc .and.		(	abs(		M_W(n,m)		) < acc	)
				end if
				!
				if(.not. test_gauge_with_acc	)	exit outer_loop
			end do
		end do outer_loop
		!
		return
	end function

!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_blas_matmul()
		logical					::	real_test, imag_test
		!
		real_test			= test_d_blas_matmul()
		imag_test			= test_z_blas_matmul()
		!
		if(.not. real_test)	call push_to_outFile("[test_blas_matmul]: failed real matrix ")
		if(.not. imag_test)	call push_to_outFile("[test_blas_matmul]: failed imag matrix ")
		!
		test_blas_matmul	= real_test .and. imag_test
		return
	end function
	!
	logical function test_d_blas_matmul()
		real(dp),	allocatable			::	A(:,:), B(:,:), C_intern(:,:), C_blas(:,:)
		real(dp)						::	acc
		character(len=130)				::	msg
		!
		allocate(	A(			smpl_size	, smpl_size		)		)
		allocate(	B(			smpl_size	, smpl_size		)		)
		allocate(	C_intern(	smpl_size	, smpl_size		)		)			
		allocate(	C_blas(		smpl_size	, smpl_size		)		)
		!
		call random_matrix(A)
		call random_matrix(B)
		!
		C_intern	=	matmul(A,B)
		call	blas_matmul(A,B,C_blas)
		!
		test_d_blas_matmul	=	.false.
		acc					= 	fp_acc
		do while ((.not. test_d_blas_matmul) .and. (acc < 1e-5_dp))
			test_d_blas_matmul	=	 is_equal_mat(acc, C_intern, C_blas)
			acc	= acc * 10.0_dp
		end do
		if(test_d_blas_matmul )	then	
			write(msg,*)	"d_blas_matmul has accuracy ",acc
			call push_to_outFile(msg)
		end if
		!
		return
	end function
	!
	logical function test_z_blas_matmul()
		complex(dp),	allocatable		::	A(:,:), B(:,:), C_intern(:,:), C_blas(:,:)
		real(dp)						::	acc
		character(len=130)				::	msg
		!
		allocate(	A(			smpl_size	, smpl_size		)		)
		allocate(	B(			smpl_size	, smpl_size		)		)
		allocate(	C_intern(	smpl_size	, smpl_size		)		)			
		allocate(	C_blas(		smpl_size	, smpl_size		)		)
		!
		call random_matrix(A)
		call random_matrix(B)
		!
		C_intern	=	matmul(A,B)
		call	blas_matmul(A,B, C_blas)
		!
		!
		test_z_blas_matmul	= .false.
		acc			= 	fp_acc
		do while ((.not. test_z_blas_matmul) .and. (acc < 1e-5_dp))
			test_z_blas_matmul	=	 is_equal_mat(acc,C_intern, C_blas)
			acc	= acc * 10_dp
		end do
		if( test_z_blas_matmul)	then
			write(msg,*)	"z_blas_matmul has accuracy ",acc
			call push_to_outFile(msg)
		end if
		!
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
		integer							::	m,n, succ_cnt
		real(dp)						::	scale, rand1, rand2
		logical							::	identity_test, wikipedia_comm, wikipedia_non_comm
		!
		!----------------------------------------------------------------------------------------------------------
		!	IDENTITY_TEST:	TEST IF RANDOM MATRIX M COMMUTES WITH IDENTITY
		!
		test_mat_comm	=	.false.
		scale			=	1023.34514_dp
		!
		allocate(	I_mat(		smpl_size, smpl_size	)		)		
		allocate(	M_mat(		smpl_size, smpl_size	)		)
		allocate(	cplx_comm(	smpl_size, smpl_size	)		)
		I_mat(:,:)	=	cmplx(0.0_dp, 0.0_dp, dp)
		M_mat(:,:)	=	cmplx(0.0_dp, 0.0_dp, dp)
		!
		!FILL MATRICES
		do n = 1, smpl_size
			do m = 1, smpl_size
				if(m==n)	I_mat(m,n)	=	cmplx(1.0,0.0_dp, dp)
				call random_number(	rand1	)			
				call random_number(	rand2	)
				M_mat(m,n)	=	cmplx( 		scale*(	rand1-.5_dp)	,	scale*(	rand2-.5_dp),	dp)
			end do
		end do
		!
		!GET COMMUTATOR
		cplx_comm	=	 matrix_comm(I_mat, M_mat)
		!
		!TEST IF COMMUTATOR IS ZERO
		succ_cnt = 0
		do n = 1, smpl_size
			do m = 1, smpl_size
				if(	 abs(	cplx_comm(m,n)	) < fp_acc		)	succ_cnt = succ_cnt + 1
			end do
		end do
		identity_test	=	( succ_cnt == smpl_size**2 )	
		!----------------------------------------------------------------------------------------------------------
		!	wikipedia_comm:	COMMUTING MATRICES TEST
		!
		A(:,:)	=	0.0_dp
		A(1,2)	=	1.0_dp
		B(:,:)	=	0.0_dp
		B(1,1)	=	1.0_dp
		B(2,2)	=	1.0_dp
		!
		real_comm		=	matrix_comm(A,B)
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
		real_comm			=	matrix_comm(A,B)
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



end program            