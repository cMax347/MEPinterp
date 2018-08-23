program test_MEPinterp
	use test_matrix_math,	only:			matrix_math_test
	use test_file_io,		only:			file_io_test

	use test_log,			only:			init_outFile, 			&
											push_to_outFile
	implicit none

	logical						::			print_also_cli,			&
											mm_success,				&
											io_success,				&
											tot_success


	integer						::			mm_array_size,				&
											io_array_size,				&
											mm_succ_cnt, mm_nTests,		&
											io_succ_cnt, io_nTests,		&
											tot_succ_cnt, tot_nTests

	character(len=120)			::			final_msg
			



	!allocation size of test arrays, matrices, etc.
	mm_array_size	= 1000
	io_array_size	= 10000




	!init the log	(with/without output to cli)
	print_also_cli	= .true.
	call init_outFile(print_also_cli, "test.log")

	tot_succ_cnt	=	0	
	tot_nTests		=	0	


	!	   * PERFORME TESTS *
	!		 --------------
	!----------------------------------------------------------------
	!				MATRIX MATH TEST
	!	
	mm_success		=	matrix_math_test(mm_array_size,		mm_succ_cnt, mm_nTests)
	tot_succ_cnt	=	tot_succ_cnt	+	 mm_succ_cnt
	tot_nTests		=	tot_nTests		+	 mm_nTests
	if(.not. mm_success)	call push_to_outFile('[test_MEPinterp]: not all matrix_math tests passed')
	call push_to_outFile("------------------------------------------------")
	!----------------------------------------------------------------
	!				FILE IO TESTS
	!
	io_success		=	file_io_test(io_array_size,		 io_succ_cnt,	io_nTests)
	tot_succ_cnt	=	tot_succ_cnt	+ 	io_succ_cnt
	tot_nTests		=	tot_nTests		+ 	io_nTests
	if(.not. io_success)	call push_to_outFile('[test_MEPinterp]: not all file_io tests passed')
	call push_to_outFile("------------------------------------------------")
	!----------------------------------------------------------------
	!
	!		PRINT FINAL MESSAGE
	!	
	call push_to_outFile("------------------------------------------------")
	call push_to_outFile("------------------------------------------------")
	call push_to_outFile("------------------------------------------------")
	write(final_msg,'(a,i4,a,i4,a)')	"[total]: passed ",tot_succ_cnt," of ",tot_nTests," tests"
	call push_to_outFile(	trim(final_msg)		)

	!
	!
	tot_success	=	mm_success .and. io_success
	if(		tot_success		)	then
		call exit(0)
	else
		call exit(0)
	end if

	!stop
end program



