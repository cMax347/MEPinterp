module test_matrix_math 
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
	logical function matrix_math_test()
		logical,			dimension(7)		::	passed
		character(len=80), dimension(7)			::	label
		character(len=10)						:: 	result
		character(len=121)						::	succ_message
		integer									::	nTests, test
		!
		nTests	= size(passed)
		!
		label(1)			=	"Levi Civitia operartor"
		passed(1)			=	test_levi_civita()
		label(2)			=	"eigen Solver"
		passed(2)			=	test_eig_solver()
		label(3)			=	"Gauge rotation"
		passed(3)			=	test_gauge_trafo()							
		label(4)			=	"Cross product"
		passed(4)			=	test_crossp()
		label(5)			=	"equality of 3D vectors"
		passed(5)			=	test_vect_equal()
		label(6)			=	"tensor to vector conversion"
		passed(6)			=	test_tens_vect()
		label(7)			=	"matrix commutator"
		passed(7)			=	test_mat_comm()
		!
		matrix_math_test	= .true.
		do test = 1, nTests
			matrix_math_test	= matrix_math_test .and. passed(test)
			if( passed(test) )	then 
				result	=	"passed"
			else
				result	=	"not passed"
			end if
			write(succ_message,'(a,a,a,a,a)')	'[test_matrix_math]: ',result,' test of "',trim(label(test)),'"'
			write(*,*)				succ_message
			!
			!WRITE TO LOG FILE
			call push_to_outFile(succ_message)
		end do
		call push_to_outFile(" ")
		!
		return
	end function







!private:------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_levi_civita()
		test_levi_civita	= 	.false.


		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_eig_solver()
		test_eig_solver		=	.false.

		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_gauge_trafo()
		test_gauge_trafo	=	.false.

		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_crossp()
		test_crossp			=	.false.

		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------
!
!
	logical function test_vect_equal()
		test_vect_equal		=	.false.

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
	logical function test_mat_comm()
		test_mat_comm		=	.false.


		return
	end function
!--------------------------------------------------------------------------------------------------------------------------------














































end module               