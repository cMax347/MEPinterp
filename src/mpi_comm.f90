module mpi_comm

#ifdef USE_MPI
	use mpi
#endif
	use constants, 		only:		dp

	implicit none


	private
	public		::		&	!parameters		
							mpi_id, mpi_root_id, mpi_nProcs, ierr,	&
							!routines
							mpi_ki_selector,						&
							!mpi_init_glob,							&
							mpi_bcast_tens,							&
							mpi_reduce_sum,						&
							mpi_allreduce_sum




!	interface mpi_init_glob
!		module procedure d2_allo_init_glob
!		module procedure d3_allo_init_glob
!		module procedure z2_allo_init_glob
!	end interface mpi_init_glob

	interface mpi_bcast_tens
		module procedure i0_mpi_bcast_tens
	end interface

	interface mpi_reduce_sum
		module procedure i0_mpi_reduce_sum
		module procedure d1_mpi_reduce_sum
		module procedure d2_mpi_reduce_sum
		module procedure d3_mpi_reduce_sum
		module procedure z2_mpi_reduce_sum
	end interface mpi_reduce_sum

	interface mpi_allreduce_sum
		module procedure i0_mpi_allreduce_sum
	end interface mpi_allreduce_sum
	



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






!private	interface subs
	!

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!
!				mpi_bcast_tens INTERFACE
!
!			performs broadcast if more then single mpi thread available
!------------------------------------------------------------------------------------------------------------

	subroutine i0_mpi_allreduce_sum( i0_loc, i0_glob)	
		integer,					intent(in)		::	i0_loc
		integer,					intent(out)		::	i0_glob
		!
		if(mpi_nProcs > 1 ) then
			call MPI_ALLREDUCE(	i0_loc, 	i0_glob,		1,		MPI_INTEGER,	MPI_SUM,	MPI_COMM_WORLD,		ierr)
			if( ierr /= 0)	stop  "[mpi_comm/i0_mpi_allreduce_sum]: (MPI_ALLREDUCE) failed"
		else
			i0_glob	=	i0_loc
		end if
		!
		return
	end subroutine

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!
!				mpi_bcast_tens INTERFACE
!
!			performs broadcast if more then single mpi thread available
!------------------------------------------------------------------------------------------------------------

	subroutine i0_mpi_bcast_tens(	i0)	
		integer,					intent(inout)		::	i0
		!
		if(mpi_nProcs > 1 ) then
			call MPI_BCAST(	i0,		1,		MPI_INTEGER,	mpi_root_id,	MPI_COMM_WORLD,		ierr			)
			if( ierr /= 0)	stop  "[mpi_comm/i0_mpi_bcast_tens]: (MPI_BCAST) failed"
		end if
		!
		return
	end subroutine

!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!
!				mpi_reduce_sum INTERFACE
!
!			performs reduction if source is allocated
!------------------------------------------------------------------------------------------------------------
	subroutine i0_mpi_reduce_sum(	loc_int,	glob_int)
		integer,							intent(in)		::	loc_int
		integer,							intent(out)		::	glob_int		
		!
		if(mpi_nProcs > 1) then
			call MPI_REDUCE(	loc_int,	glob_int,	1,	MPI_INTEGER,	MPI_SUM,	mpi_root_id,	MPI_COMM_WORLD, ierr)
			if(ierr /= 0)	stop "[mpi_comm/i0_mpi_reduce_sum]: (MPI_REDUCE) failed"
		else
			glob_int = loc_int
		end if
		!
		return
	end subroutine



	subroutine d1_mpi_reduce_sum(	loc_tens, glob_tens)
		real(dp),							intent(in)		::	loc_tens(:)
		real(dp),		allocatable,		intent(inout)	::	glob_tens(:)
		integer												::	package_size
		!	
		allocate(	glob_tens(	size(loc_tens,1))		)
		write(*,*)	"[d1_mpi_reduce_sum]: allocated 1d real(dp) tensor with size ",size(loc_tens,1)
		!
		!
		if(	mpi_nProcs >1) then
			package_size	=	size(loc_tens)
			call MPI_REDUCE(	loc_tens,	glob_tens,	package_size,	MPI_DOUBLE_PRECISION,	MPI_SUM, mpi_root_id, MPI_COMM_WORLD, ierr)
			if(ierr /= 0)	stop "[mpi_comm/d1_mpi_reduce_sum]: (MPI_REDUCE) failed"
		else
			glob_tens	=	loc_tens
		end if
		!
		return
	end subroutine


	subroutine d2_mpi_reduce_sum( loc_tens, glob_tens)	
		real(dp),							intent(in)		::	loc_tens(:,:)
		real(dp),		allocatable,		intent(inout)	::	glob_tens(:,:)
		integer												::	package_size
		!
		allocate(		glob_tens(	size(loc_tens,1), size(loc_tens,2)	))
		glob_tens=	0.0_dp
		!
		!
		if(	mpi_nProcs > 1) then
			package_size	=	size(loc_tens,1)	* size(loc_tens,2)
			call MPI_REDUCE(	loc_tens,  glob_tens, package_size, 	MPI_DOUBLE_PRECISION, MPI_SUM, mpi_root_id, MPI_COMM_WORLD,	ierr)
			if(ierr /= 0)	stop "[mpi_comm/d2_mpi_reduce_sum]: (MPI_REDUCE) failed"
		else
			glob_tens	=	loc_tens
		end if
		!
		!
		return
	end subroutine	


	subroutine d3_mpi_reduce_sum( loc_tens, glob_tens)	
		real(dp),							intent(in)		::	loc_tens(:,:,:)
		real(dp),		allocatable,		intent(inout)	::	glob_tens(:,:,:)
		integer												::	package_size
		!
		allocate(		glob_tens(	size(loc_tens,1), size(loc_tens,2), size(loc_tens,3)	))
		glob_tens	= 0.0_dp
		!
		!
		if( mpi_nProcs > 1) then
			package_size	=	size(loc_tens,1)	* size(loc_tens,2) * size(loc_tens,3)
			call MPI_REDUCE(	loc_tens,  glob_tens, package_size, 	MPI_DOUBLE_PRECISION, MPI_SUM, mpi_root_id, MPI_COMM_WORLD,	ierr)
			if(ierr /= 0)	stop "[mpi_comm/d3_mpi_reduce_sum]: (MPI_REDUCE) failed"
		else
			glob_tens	=	loc_tens
		end if
		!
		!
		return
	end subroutine	


	subroutine z2_mpi_reduce_sum( loc_tens, glob_tens)	
		complex(dp),							intent(in)		::	loc_tens(:,:)
		complex(dp),		allocatable,		intent(inout)	::	glob_tens(:,:)
		integer													::	package_size
		!
		allocate(		glob_tens(	size(loc_tens,1), size(loc_tens,2)	))
		glob_tens	=	cmplx(	0.0_dp, 0.0_dp,	dp)
		!
		!
		if( mpi_nProcs > 1) then
			package_size	=	size(loc_tens,1)	* size(loc_tens,2)
			call MPI_REDUCE(	loc_tens,  glob_tens, package_size, 	MPI_DOUBLE_COMPLEX , MPI_SUM, mpi_root_id, MPI_COMM_WORLD,	ierr)
		if(ierr /= 0)	stop "[mpi_comm/z2_mpi_reduce_sum]: (MPI_REDUCE) failed"
		else
			glob_tens	=	loc_tens
		end if
		!
		!
		return
	end subroutine	











end module mpi_comm