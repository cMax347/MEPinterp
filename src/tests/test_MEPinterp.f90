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


	integer						::			mm_succ_cnt, mm_nTests,		&
											io_succ_cnt, io_nTests,		&
											tot_succ_cnt, tot_nTests

	character(len=120)			::			final_msg
			
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
	mm_success		=	matrix_math_test(mm_succ_cnt, mm_nTests)
	tot_succ_cnt	=	tot_succ_cnt	+	 mm_succ_cnt
	tot_nTests		=	tot_nTests		+	 mm_nTests
	!----------------------------------------------------------------
	!				FILE IO TESTS
	!
	io_success		=	file_io_test(io_succ_cnt,	io_nTests)
	tot_succ_cnt	=	tot_succ_cnt	+ 	io_succ_cnt
	tot_nTests		=	tot_nTests		+ 	io_nTests
	!
	
	!----------------------------------------------------------------
	!
	!		PRINT FINAL MESSAGE
	!	
	call push_to_outFile("------------------------------------------------")
	write(final_msg,'(a,i4,a,i4)')	"[total]: passed ",tot_succ_cnt," of ",tot_nTests," tests"
	call push_to_outFile(	trim(final_msg)		)

	!
	!
	tot_success	=	mm_success .and. io_success
	if( .not. tot_success	)	stop "not all tests where passsed"


	!stop
end program



