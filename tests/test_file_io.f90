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

	integer					::		smpl_size

	contains


!public:
	logical function file_io_test(test_arr_size, succ_cnt, nTests)
		integer,			intent(in)			::	test_arr_size
		integer,			intent(out)			::	succ_cnt, nTests
		logical,			allocatable			::	passed(:)
		character(len=80), 	allocatable			::	label(:)
		!
		smpl_size		=	test_arr_size
		file_io_test	= 	.false.
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
		!
		!		RETURN SUCCESS
		!
		file_io_test	=	(succ_cnt == nTests)
		return
	end function










!private
	logical function test_bin_en_io()
		real(dp),		allocatable			::	en_orig(:), en_clone(:)
		real(dp)							::	rand
		integer								::	qi_idx, idx


		allocate(	en_orig(smpl_size)		)		
		allocate(	en_clone(smpl_size)		)
		!
		!	try to write at random index								
		call random_number(rand)
		qi_idx	= 100 + int(	rand * 100	)
		!	a random array
		call d_rand_arr(	en_orig		)
		!
		!											WRITE: orignal 			READ: clone
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
		!
		return
	end function

	subroutine d_rand_arr(v)
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