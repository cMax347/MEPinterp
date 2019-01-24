program MEPinterp
!
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
#ifdef USE_MPI	
	use mpi
#endif
	use mpi_community,			only:		mpi_root_id, mpi_id, mpi_nProcs, ierr
	use input_paras,			only:		init_parameters,	&
											plot_bands
	use core,					only:		core_worker
	use band_calc,				only:		band_worker
	implicit none
	!
	!
	!MPI INIT
	mpi_root_id = 	0
	mpi_id		=	0
	mpi_nProcs	=	1

#ifdef USE_MPI
	call MPI_INIT( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, 	mpi_id			, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, 	mpi_nProcs		, ierr)
#endif
	write(*,'(a,i3,a,a,a)')	'[#',mpi_id,';main/',cTIME(time()),']:	welcome to mepInterp'
	!
	!
    !BODY
    if( init_parameters()	) then
		if(	plot_bands	) then
			call band_worker()
		else
			call core_worker()			
		end if
	else
		write(*,'(a,i3,a)')		'[#',mpi_id,';main]: input file not found, by'
	end if
	!
	!
	!FINALIZE
#ifdef USE_MPI	
	call MPI_FINALIZE(ierr)
#endif
	write(*,'(a,i3,a,a,a)')	'[#',mpi_id,';main/',cTIME(time()),']:	all done, by by'	
	!
	!
	stop
end program 