module test_file_io
	use file_io,		only:						mpi_read_k_mesh, mpi_read_tb_basis,	& 
													write_en_binary, read_en_binary,	&
													write_en_global,					&
													write_mep_tensors
    use test_log,		only:				push_to_outFile                         


	implicit none

	private
	public					::		test_file_io

	contains



	logical function test_file_io()

		test_file_io	= .false.

		call push_to_outFile("[test_file_io]: tests not implemented yet")


		return
	end function

end module