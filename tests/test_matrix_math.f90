module test_matrix_math 
	use parameters,		only:				dp, fp_acc
	use matrix_math,	only:				my_Levi_Civita,                         &
                                            zheevr_wrapper, zheevd_wrapper,         & 
                                            uni_gauge_trafo,                        &
                                            is_equal_vect,                          &
                                            convert_tens_to_vect,                   &
                                            blas_matmul,							&
                                            matrix_comm             

    use test_log,		only:				push_to_outFile, write_test_results                         

	implicit none


	private
	public				::				matrix_math_test


	interface random_matrix
		module procedure ds_rand_mat
		module procedure zh_rand_mat	
	end interface random_matrix

	interface is_equal_mat
		module procedure d_is_equal_mat
		module procedure z_is_equal_mat
	end interface is_equal_mat


	integer				::				smpl_size

	contains

!public:
	logical function matrix_math_test(test_mat_size, succ_cnt, nTests)
		integer,			intent(in)		::	test_mat_size
		integer,			intent(out)		::	succ_cnt, nTests
		logical,		   allocatable		::	passed(:)
		character(len=80), allocatable		::	label(:)
		!
		smpl_size	= test_mat_size
		!
		matrix_math_test	= .false.
		nTests	= 5
		succ_cnt= 0
		allocate(	passed(nTests)	)		
		allocate(	label(nTests)	)		
		!
		!
		label(		1	)			=	"Levi Civitia operartor"
		passed(		1	)			=	test_levi_civita()
		!-----------------------------------------------------------------
		label(		2	)			=	"equality of  vectors"
		passed(		2	)			=	test_vect_equal()
		!-----------------------------------------------------------------
		label(		3	)			=	"Gauge rotation"
		passed(		3	)			=	test_gauge_trafo()				
		!-----------------------------------------------------------------			
		label(		4	)			=	"matmul with Blas"
		passed(		4	)			=	test_blas_matmul()
		!-----------------------------------------------------------------			
		label(		5	)			=	"matrix commutator"
		passed(		5	)			=	test_mat_comm()
		!
		!
		!	WRITE RESULTS TO LOG
		!
		call write_test_results(	'matrix_math_test', passed, label, succ_cnt, nTests	)
		!
		matrix_math_test	=	( succ_cnt	==	nTests)
		!
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
		real(dp), 	allocatable				::	a(:), b(:), a_cpy(:)
		integer								::	smpl_size
		!
		smpl_size	= 1000
		allocate( a(		smpl_size)	)
		allocate( b(		smpl_size)	)		
		allocate( a_cpy(	smpl_size)	)
		!
		call d_rand_vec(a)
		call d_rand_vec(b)
		b(1)		=	a(1)-1e-13_dp
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
	logical function test_tens_vect()
		test_tens_vect		=	.false.

		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_gauge_trafo()
		complex(dp),	allocatable			::	M_W(:,:),				&
												U(:,:)
		real(dp),		allocatable			::	eigVal(:)
		integer								::	n,m
		!
		test_gauge_trafo	=	.false.
		!
		allocate(	M_W(	smpl_size, smpl_size	)		)
		allocate(	U(		smpl_size, smpl_size	)		)
		allocate(	eigVal(	smpl_size		)		)
		!
		!random hermitian matrix
		call zh_rand_mat(M_W)
		!GET EIGVAL & EIGVEC
		U	=	M_W
		call zheevd_wrapper(U,	eigVal)
		write(*,*)	"[test_gauge_trafo]: solved random hermitian matrix with eigvals:"
		write(*,*)	eigVal
		!
		!rotate M_W with U should yield diagonal matrix (matrix in H gauge)
		call uni_gauge_trafo(U, M_W)
		!
		!
		!todo: check now if M_H is diagonal with eival on diago
		test_gauge_trafo	= .true.
		outer_loop:	do m = 1, size(M_W,2)
			do n = 1, size(M_W,1)
			
				if(	n==m	) then
					test_gauge_trafo =	test_gauge_trafo .and. 		(	abs(	M_W(n,n)-eigVal(n) 	) < fp_acc	)		
				else
					test_gauge_trafo =	test_gauge_trafo .and.		(	abs(		M_W(n,m)		) < fp_acc	)
				end if
				!
				!
				if(.not. test_gauge_trafo	)	then
					write(*,*)	"[test_gauge_trafo]:	unexpected value at  |M(n=",n,',m=',m,')|= ', abs(M_W(m,n))
					if(n==m)	write(*,*)	"[test_gauge_trafo]:	eigen value=",eigVal(n)
					!exit outer_loop
				end if
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
		C_blas		=	blas_matmul(A,B)
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
		C_blas		=	blas_matmul(A,B)
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








!helpers:
!------------------------------------------------------------------------------------------------------
	subroutine ds_rand_mat(Mat)
		real(dp), intent(out)	::	Mat(:,:)
		integer						::	n, m
		real(dp)					::	re_rand
		!
		!
		do n = 1, size(Mat,2)
			do m = 1, size(Mat,1)
				call random_number(re_rand)	
				!
				re_rand	= (re_rand - .5_dp) 							
				!
				Mat(n,m)	=	re_rand
				if(n<m) then
					Mat(m,n)	=	re_rand
				end if
				!
				!
			end do
		end do
		!
		return
	end subroutine	


	subroutine zh_rand_mat(Mat)
		complex(dp), intent(out)	::	Mat(:,:)
		integer						::	n, m
		complex(dp)					::	rand
		real(dp)					::	re_rand, im_rand
		!
		do n = 1, size(Mat,2)
			do m = 1, size(Mat,1)
				call random_number(re_rand)
				re_rand	= (re_rand - .5_dp) 		
				!
				if(	n==m ) then
					Mat(n,n)	=	re_rand
				else
					call random_number(im_rand)	
					im_rand	= (im_rand - .5_dp)					
					rand	= cmplx(	re_rand, im_rand,	dp)
					!
					Mat(n,m)	=	rand
					Mat(m,n)	= conjg(rand)
				end if
				!
				!
			end do
		end do
		!
		return
	end subroutine	
!------------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------------
	logical function	d_is_equal_mat(acc, A, B)
		real(dp),	intent(in)			::	acc
		real(dp),	intent(in)			::	A(:,:), B(:,:)
		real(dp)						::	delta
		integer							::	n, m
		!
		d_is_equal_mat	=	(	size(A,1) == size(B,1) .and. size(A,2) == size(B,2) )
		!
		if( d_is_equal_mat	) then

			loop_out: do n = 1, size(A,2)
				do m = 1, size(A,1)
					delta	=	abs(	A(m,n) - B(m,n)	)
					!
					d_is_equal_mat	= d_is_equal_mat .and. (delta < acc)
					!
					if(.not. d_is_equal_mat) exit loop_out
				end do
			end do loop_out
		end if
		!
		return
	end function

	logical function	z_is_equal_mat(acc, A, B)
		real(dp),		intent(in)		::	acc
		complex(dp),	intent(in)		::	A(:,:), B(:,:)
		real(dp)						::	delta
		integer							::	n, m
		character(len=120)				::	msg
		!
		z_is_equal_mat	=	(	size(A,1) == size(B,1) .and. size(A,2) == size(B,2) )
		!
		if( z_is_equal_mat	) then

			loop_out: do n = 1, size(A,2)
				loop_in: do m = 1, size(A,1)
					delta	=	abs(	A(m,n) - B(m,n)	)
					!
					z_is_equal_mat	= z_is_equal_mat .and. (delta < acc)
					!
					if(.not. z_is_equal_mat)then 
						!write(msg,'(a,f4.2,a,f4.2,a,f4.2,a,f4.2,a)')"[z_is_equal_mat]: (",real(A(m,n)),"+i",imag(A(m,n)) ,") vs (", real(B(m,n)),"+i",imag(B(m,n)),")"
						call push_to_outFile(msg)
						exit loop_out
					end if
				end do loop_in
			end do loop_out
		else
			write(msg,*)	"[z_is_equal_mat]: sizes to not compare"
			call push_to_outFile(msg)
		end if
		!
		return
	end function
!------------------------------------------------------------------------------------------------------

	subroutine d_rand_vec(v)
		real(dp),	intent(out)		::	v(:)
		integer						:: idx
		real(dp)					:: rand, scal
		!
		call random_number(scal)
		scal	=	scal * size(v,1)
		!
		do idx = 1, size(v,1)
			call random_number(rand)
			rand = (rand-.5_dp)	* scal
			!
			!
			v(idx)	= rand
		end do
		!
		return
	end subroutine





































end module               