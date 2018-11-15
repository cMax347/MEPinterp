program test_file_io
	use constants, 		only:						dp,	fp_acc	
	use file_io,		only:						read_tb_basis,					& 
													write_en_binary, read_en_binary,	&
													write_en_global
    use helpers,		only:						my_exit,							&
    												write_test_results, push_to_outFile,&
    												random_vector                        


	implicit none

	

	integer									::		smpl_size,		& 
													nTests, succ_cnt
	logical									::		all_passed	
	logical,			allocatable			::		passed(:)
	character(len=80), 	allocatable			::		label(:)
	!
	smpl_size		=	1000
	all_passed		= 	.false.
	nTests			= 	1
	succ_cnt		= 	0
	allocate(	label(nTests)	)
	allocate(	passed(nTests)	)		
	!
	!		GIVE EACH TEST A LABEL THEN PERFORM TEST				
	!
	label(		1	)			=	"binary energy I/O"
	passed(		1	)			=	test_bin_en_io()
	!-----------------------------------------------------------------
	!
	!		WRITE LOG FILE
	call write_test_results("file_io_test", passed, label, succ_cnt, nTests)
	call push_to_outFile("------------------------------------------------------")
	call push_to_outFile("")
	!
	!		RETURN SUCCESS
	!
	all_passed	=	(succ_cnt == nTests)
	call my_exit(all_passed)




contains


	!private
	logical function test_bin_en_io()
		real(dp),		allocatable			::	en_orig(:), en_clone(:)
		real(dp)							::	rand
		integer								::	qi_idx, idx
		!
		allocate(	en_orig(smpl_size)		)		
		allocate(	en_clone(smpl_size)		)
		!							
		call random_number(rand)
		qi_idx	= 100 + int(	rand * 100	)
		!
		call random_vector(	en_orig		)					
		call write_en_binary(	qi_idx, 	en_orig		)
		call read_en_binary(	qi_idx,		en_clone	)
		!
		!	COMPARE ORIG & CLONE
		!
		test_bin_en_io	=	.true.
		do idx = 1, smpl_size
			test_bin_en_io	=	test_bin_en_io .and. 	(	abs( en_orig(idx) - en_clone(idx) )		<	fp_acc	)
		end do
		!
		return
	end function




end program















