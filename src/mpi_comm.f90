module mpi_comm

#ifdef USE_MPI
	use mpi
#endif

	implicit none


	private
	public		::		&	!parameters		
							mpi_id, mpi_root_id, mpi_nProcs, ierr,	&
							!routines
							mpi_ki_selector



	save
	!MPI
	integer						::	mpi_id, mpi_root_id, mpi_nProcs, ierr

contains



!public
	



	logical pure function mpi_ki_selector(ki_request, num_kpts)
		!		map the kpts onto the mpi-threads running
		integer,		intent(in)		::		ki_request, num_kpts 
		integer							::		ki_todo
		!
		mpi_ki_selector = .false.
		!
		if(		mpi_nProcs == 1		) then
			mpi_ki_selector	=	.true.
		else
			loop_todos: do ki_todo = mpi_id +1, num_kpts, mpi_nProcs
				if(	ki_request == ki_todo) then
					mpi_ki_selector = .true.
					exit loop_todos
				end if
			end do loop_todos
		end if
		!
		return
	end function

	!-----------------------------------------------------------------------------------
	!-----------------------------------------------------------------------------------



!private






end module mpi_comm