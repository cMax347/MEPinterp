program MEPinterp
	use mpi
	use parameters,				only:		init_parameters,									&
											dp,													&
											mpi_root_id, mpi_id, mpi_nProcs, ierr
	use mep_niu,				only:		mep_interp
	use file_io,				only:		write_mep_tensor
	implicit none


	real(dp),		dimension(3,3)	::		mep_tens_loc, 										&
											mep_tens_glob


	!MPI INIT
	mpi_root_id = 0
	call MPI_INIT( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, 	mpi_id			, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, 	mpi_nProcs		, ierr)

    !read the input file
    call init_parameters()

	!do the interpolation
	call mep_interp(mep_tens_loc)			


	!
	call MPI_REDUCE(	mep_tens_loc, mep_tens_glob,	 9,	MPI_DOUBLE_PRECISION,	MPI_SUM	, 	mpi_root_id,	MPI_COMM_WORLD, ierr)
	if(mpi_id == mpi_root_id) 	call write_mep_tensor(mep_tens_glob)


	call MPI_FINALIZE(ierr)
	stop

end program 