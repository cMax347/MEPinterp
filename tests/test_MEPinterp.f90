program test_MEPinterp
	use test_matrix_math,	only:			matrix_math_test
	use test_file_io,		only:			file_io_test

	use test_log,			only:			init_outFile
	implicit none

	logical						::			matrix_math_success
			
	!init the log
	call init_outFile("test.log")

	!performe tests
	matrix_math_success		=	matrix_math_test()

	



	stop
end program



