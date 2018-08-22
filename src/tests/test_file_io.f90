module test_file_io
	use parameters,		only:						dp,	fp_acc							
	use file_io,		only:						mpi_read_k_mesh, mpi_read_tb_basis,	& 
													write_en_binary, read_en_binary,	&
													write_en_global,					&
													write_mep_tensors
    use test_log,		only:						write_test_results                        


	implicit none

	private
	public					::		file_io_test

	contains


!public:
	logical function file_io_test(succ_cnt, nTests)
		integer,			intent(out)			::	succ_cnt, nTests
		logical,			dimension(1)		::	passed
		character(len=80), 	dimension(1)		::	label
		!
		!		GIVE EACH TEST A LABEL THEN PERFORM TEST				
		!
		label(		1	)			=	"binary energy I/O"
		passed(		1	)			=	test_bin_en_io()
		!-----------------------------------------------------------------
		!
		!		WRITE LOG FILE
		call write_test_results("file_io_test", passed, label, succ_cnt, nTests)
		!
		!		RETURN SUCCESS
		!
		file_io_test	=	(succ_cnt == nTests)
		return
	end function










!private
	logical function test_bin_en_io()
		real(dp),		allocatable			::	en_orig(:), en_clone(:)
		real(dp)							::	rand1, rand2
		integer								::	size, qi_idx, idx


		size = 20
		allocate(	en_orig(size)		)		
		allocate(	en_clone(size)		)
		!
		!	FILL original	(RANDOM VALUES)
		!
		call random_number(rand1)
		qi_idx	= 100 + int(	rand1 * 100	)
		do idx	= 1, size
			call random_number(rand1)
			call random_number(rand2)
			en_orig(idx)	=	 (	rand1-.5_dp	) 	* rand2 * size**2
		end do
		!
		!	WRITE: orignal 			READ: clone
		!
		write(*,*)	"en_orig=",en_orig
		call write_en_binary(	qi_idx, 	en_orig(1:size)		)
		call read_en_binary(	qi_idx,		en_clone(1:size)	)
		!
		!	COMPARE ORIG & CLONE
		!
		test_bin_en_io	=	.true.
		do idx = 1, size
			test_bin_en_io	=	test_bin_en_io .and. 	(	abs( en_orig(idx) - en_clone(idx) )		<	fp_acc	)
		end do
		!
		!
		return
	end function




end module