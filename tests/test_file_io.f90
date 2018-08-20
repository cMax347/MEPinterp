module test_file_io
!	use file_io,		only:						mpi_read_k_mesh, mpi_read_tb_basis,	& 
!													write_en_binary, read_en_binary,	&
!													write_en_global,					&
!													write_mep_tensors
    use test_log,		only:				push_to_outFile                         


	implicit none

	private
	public					::		file_io_test

	contains



	logical function file_io_test()

		file_io_test	= .false.

		call push_to_outFile("[file_io_test]: tests not implemented yet")


		return
	end function

end module