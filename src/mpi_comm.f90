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
							mpi_reduce_tens




!	interface mpi_init_glob
!		module procedure d2_allo_init_glob
!		module procedure d3_allo_init_glob
!		module procedure z2_allo_init_glob
!	end interface mpi_init_glob

	interface mpi_reduce_tens
		module procedure i0_mpi_reduce_tens
		module procedure d2_mpi_reduce_tens
		module procedure d3_mpi_reduce_tens
		module procedure z2_mpi_reduce_tens
	end interface mpi_reduce_tens






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
!				mpi_reduce_tens INTERFACE
!
!			performs reduction if source is allocated
!------------------------------------------------------------------------------------------------------------
	subroutine i0_mpi_reduce_tens(	loc_int,	glob_int)
		integer,							intent(in)		::	loc_int
		integer,							intent(out)		::	glob_int		
		!
		if(mpi_nProcs > 1) then
			call MPI_REDUCE(	loc_int,	glob_int,	1,	MPI_INTEGER,	MPI_SUM,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		else
			glob_int = loc_int
		end if
		!
		return
	end subroutine

	subroutine d2_mpi_reduce_tens( loc_tens, glob_tens)	
		real(dp),		allocatable,		intent(inout)	::	loc_tens(:,:)
		real(dp),		allocatable,		intent(inout)	::	glob_tens(:,:)
		integer												::	package_size
		!
		if(	allocated(loc_tens)	)	then
			allocate(		glob_tens(	size(loc_tens,1), size(loc_tens,2)	))
			glob_tens=	0.0_dp

			if(	mpi_nProcs > 1) then
				package_size	=	size(loc_tens,1)	* size(loc_tens,2)
				call MPI_REDUCE(	loc_tens,  glob_tens, package_size, 	MPI_DOUBLE_PRECISION, MPI_SUM, mpi_root_id, MPI_COMM_WORLD,	ierr)
			else
				glob_tens	=	loc_tens
			end if
		end if
		!
		return
	end subroutine	


	subroutine d3_mpi_reduce_tens( loc_tens, glob_tens)	
		real(dp),		allocatable,		intent(inout)	::	loc_tens(:,:,:)
		real(dp),		allocatable,		intent(inout)	::	glob_tens(:,:,:)
		integer												::	package_size
		!
		if(	allocated(loc_tens)	)	then
			allocate(		glob_tens(	size(loc_tens,1), size(loc_tens,2), size(loc_tens,3)	))
			glob_tens	= 0.0_dp
			!
			if( mpi_nProcs > 1) then
				package_size	=	size(loc_tens,1)	* size(loc_tens,2) * size(loc_tens,3)
				call MPI_REDUCE(	loc_tens,  glob_tens, package_size, 	MPI_DOUBLE_PRECISION, MPI_SUM, mpi_root_id, MPI_COMM_WORLD,	ierr)
			else
				glob_tens	=	loc_tens
			end if
		end if
		!
		return
	end subroutine	


	subroutine z2_mpi_reduce_tens( loc_tens, glob_tens)	
		complex(dp),		allocatable,		intent(inout)	::	loc_tens(:,:)
		complex(dp),		allocatable,		intent(inout)	::	glob_tens(:,:)
		integer													::	package_size
		!
		if(	allocated(loc_tens)	)	then
			allocate(		glob_tens(	size(loc_tens,1), size(loc_tens,2)	))
			glob_tens	=	cmplx(	0.0_dp, 0.0_dp,	dp)
			!
			if( mpi_nProcs > 1) then
				package_size	=	size(loc_tens,1)	* size(loc_tens,2)
				call MPI_REDUCE(	loc_tens,  glob_tens, package_size, 	MPI_DOUBLE_COMPLEX , MPI_SUM, mpi_root_id, MPI_COMM_WORLD,	ierr)
			else
				glob_tens	=	loc_tens
			end if
		end if
		!
		return
	end subroutine	
























!
!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!
!				allo_glob_arr INTERFACE
!------------------------------------------------------------------------------------------------------------




!	subroutine	d2_allo_init_glob(	loc_arr, glob_arr)
!		!	loc:	sum over local k-pts 		(k-pts of this mpi node)
!		!	glob:	sum over complete 1st BZ	(k-pts on all mpi nodes)
!		real(dp),		allocatable,	intent(in)			::	loc_arr(:,:)
!		real(dp),		allocatable,	intent(inout)		::	glob_arr(:,:)
!		!
!		if(allocated(loc_arr))		then
!			!
!			allocate(		glob_arr(	size(loc_arr,1),	size(loc_arr,2)					))
!			if(mpi_id == mpi_root_id) 	then
!				if(	 mpi_nProcs ==	1 ) then	!SERIAL MODE - cpy local to global
!					glob_arr	=	loc_arr
!				else if(mpi_nprocs > 1) then	!PARALLEL MODE - init global containers to zero
!					glob_arr	=	0.0_dp
!				else
!					stop "[d_allo_glob_arr]:	detected negative mpi_nprocs!"
!				end if
!			end if
!		end if
!		!
!		return
!	end subroutine
!
!
!	subroutine d3_allo_init_glob(	loc_arr, glob_arr)
!		!	loc:	sum over local k-pts 		(k-pts of this mpi node)
!		!	glob:	sum over complete 1st BZ	(k-pts on all mpi nodes)
!		real(dp),		allocatable,	intent(in)			::	loc_arr(:,:,:)
!		real(dp),		allocatable,	intent(inout)		::	glob_arr(:,:,:)
!		!
!		if(allocated(loc_arr))		then
!			!
!			allocate(	glob_arr(	size(loc_arr,1),	size(loc_arr,2), 	size(loc_arr,3)		))
!			if(mpi_id == mpi_root_id) 	then
!				if(	 mpi_nProcs ==	1 ) then	!SERIAL MODE - cpy local to global
!					glob_arr	=	loc_arr
!				else if(mpi_nprocs > 1) then	!PARALLEL MODE - init global containers to zero
!					glob_arr	=	0.0_dp
!				else
!					stop "[d_allo_glob_arr]:	detected negative mpi_nprocs!"
!				end if
!			end if
!		end if
!		!
!		return
!	end subroutine
!
!
!	subroutine	z2_allo_init_glob(	loc_arr, glob_arr)
!		!	loc:	sum over local k-pts 		(k-pts of this mpi node)
!		!	glob:	sum over complete 1st BZ	(k-pts on all mpi nodes)
!		complex(dp),		allocatable,	intent(in)			::	loc_arr(:,:)
!		complex(dp),		allocatable,	intent(inout)		::	glob_arr(:,:)
!		!
!		if(allocated(loc_arr))		then
!			allocate(		glob_arr(	size(loc_arr,1),	size(loc_arr,2)						))
!			!
!			if(mpi_id == mpi_root_id) 	then
!				if(	 mpi_nProcs ==	1 ) then	!SERIAL MODE - cpy local to global
!					glob_arr	=	loc_arr
!				else if(mpi_nprocs > 1) then	!PARALLEL MODE - init global containers to zero
!					glob_arr	=	cmplx(	0.0_dp,	0.0_dp, dp)
!				else
!					stop "[d_allo_glob_arr]:	detected negative mpi_nprocs!"
!				end if
!			end if
!		end if
!		!
!		return
!	end subroutine
!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------






end module mpi_comm