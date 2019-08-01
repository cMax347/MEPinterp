program MEPinterp
!
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
	use omp_lib
#ifdef USE_MPI	
	use mpi
#endif
	use mpi_community,			only:		mpi_root_id, mpi_id, mpi_nProcs, ierr
	use input_paras,			only:		init_parameters,	&
											my_mkdir,			&
											plot_bands
	use core,					only:		core_worker
	use band_calc,				only:		band_worker
	implicit none
	!
	real 	::	T_start, T_finish
	integer	::	t_hour, t_min, t_sec
    call cpu_time(T_start)
	!
	!MPI INIT
	mpi_root_id = 	0
	mpi_id		=	0
	mpi_nProcs	=	1
	!
#ifdef USE_MPI
	call MPI_INIT( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, 	mpi_id			, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, 	mpi_nProcs		, ierr)
#endif
	write(*,'(a,i7.7,a,a,a)')	'[#',mpi_id,';main/',cTIME(time()),']:	welcome to mepInterp'
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
		write(*,'(a,i7.7,a)')		'[#',mpi_id,';main]: input file not found, by'
	end if
	!
	!
	!FINALIZE
	write(*,'(a,i7.7,a,a,a)')	'[#',mpi_id,';main/',cTIME(time()),']:	all done, by by'	
	!
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifdef USE_MPI	
	call MPI_FINALIZE(ierr)
#endif
	!
	!REPORT WALL TIME
	if(mpi_id==mpi_root_id) then
		write(*,*)	'~~~~~'
		
		call cpu_time(T_finish)
		t_hour	=	int(	(	(T_finish-T_start) 			-		mod(T_finish-T_start,60.0**2)		)	)	/ 60**2		
		t_min	=	int(	(mod(T_finish-T_start,60.0**2) 	- 	mod(mod(T_finish-T_start,60.0**2),60.0)) 	)	/ 	60		
		t_sec	=	int(	mod(mod(T_finish-T_start,60.0**2),60.0)												)
		write(*,'(a,i7.7,a,a,a,f15.3,a,i2.2,a,i2.2,a,i2.2,a)')	'[#',mpi_id,';main/',cTIME(time()),&
																']:	approximated Wall time: ', T_finish-T_start, &
																	" seconds. (hh:mm:ss::~~~",t_hour,":",t_min,":",t_sec,")"

	end if
	!
	stop
end program 