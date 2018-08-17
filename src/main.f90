program MEPinterp
	use mpi
	use parameters,				only:		init_parameters,									&
											mpi_root_id, mpi_id, mpi_nProcs, ierr
	use mep_niu,				only:		mep_interp
	implicit none


	!MPI INIT
	mpi_root_id = 0
	call MPI_INIT( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, 	mpi_id			, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, 	mpi_nProcs		, ierr)

    !read the input file
    call init_parameters()

	!do the interpolation
	call mep_interp()			


	call MPI_FINALIZE(ierr)
	stop

end program 