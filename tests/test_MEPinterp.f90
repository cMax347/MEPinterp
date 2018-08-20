program test_MEPinterp
	use test_matrix_math,	only:			test_matrix_math
	use test_file_io,		only:			test_file_io

	use test_log,			only:			init_outFile
	implicit none


			
	!init the log
	call init_outFile("test.log")

	!performe tests
	call test_matrix_math()

	



	stop
end program



