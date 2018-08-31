program MEPinterp
	use mpi
	use constants,				only:		mpi_root_id, mpi_id, mpi_nProcs, ierr
	use input_paras,			only:		init_parameters,	&
											plot_bands
	use mep_niu,				only:		mep_worker
	use band_calc,				only:		band_worker
	implicit none
	!
	!
	!MPI INIT
	mpi_root_id = 0
	call MPI_INIT( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, 	mpi_id			, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, 	mpi_nProcs		, ierr)
    !
    !read the input file
    call init_parameters()
    !
	if(	plot_bands	) then
		call band_worker()
	else
		call mep_worker()			
	end if
	!
	call MPI_FINALIZE(ierr)
	stop
end program 