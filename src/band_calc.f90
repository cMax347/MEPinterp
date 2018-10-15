module band_calc
!
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
#ifdef USE_MPI
	use mpi
#endif
	use constants,		only:			dp, mpi_root_id, mpi_id, mpi_nProcs, ierr
	use input_paras,	only:			seed_name, a_latt
	use k_space,		only:			get_recip_latt										
	use file_io,		only:			read_kptsgen_pl_file,							&
										read_tb_basis,									&
										write_en_binary, 								&
										write_en_global
	use wann_interp,	only:			get_wann_interp									
	implicit none


	private
	public		::						band_worker


contains

	subroutine band_worker()
		real(dp),		allocatable			::	rel_kpts(:,:), en_k(:), R_vect(:,:)
		real(dp)							::	recip_latt(3,3)
		integer								::	num_kpts, ki, k_per_mpi
		complex(dp),	allocatable			::	H_tb(:,:,:), r_tb(:,:,:,:), 			&
												A_ka(:,:,:), Om_kab(:,:,:,:),				&
												V_ka(:,:,:)
		logical								::	do_gauge_trafo
		!
		if( read_kptsgen_pl_file(rel_kpts)	) then
			!
			!	get k-space
			recip_latt	= get_recip_latt()
			num_kpts	= size(rel_kpts,2)
			k_per_mpi	= 0
			!
			!	get the data
			call read_tb_basis(seed_name, R_vect, H_tb, r_tb)
			call k_space_allocator(H_tb, r_tb, en_k, V_ka, A_ka, Om_kab)
			!
			!
			!	do the work
#ifdef USE_MPI
			call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
			write(*,'(a,i3,a,a,a)')	'[#',mpi_id,';band_worker/',cTIME(time()),	']:	start interpolating...'
			!
			!
			do_gauge_trafo	= .false. !eigenvalues are gauge independent
			!
			do ki = mpi_id + 1, num_kpts,	mpi_nProcs
				call get_wann_interp(do_gauge_trafo, H_tb, r_tb, a_latt, recip_latt, R_vect, rel_kpts(:,ki), 	en_k, V_ka, A_ka, Om_kab )
				call write_en_binary(ki,en_k)
				k_per_mpi	= k_per_mpi + 1
			end do
			write(*,'(a,i3,a,a,a,i10,a,i10,a)')	'[#',mpi_id,';band_worker/',cTIME(time()),				&
													']:	...finished, interpolated ',k_per_mpi,' of ',num_kpts,' kpts'
			!
			!
			!	write the results
#ifdef USE_MPI
			call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif			
			if(mpi_id == mpi_root_id)	then
				call write_en_global(rel_kpts)
				write(*,*)'---------------------------------------------------------------------------------------------'
			end if
		else
			stop 'for bandstructure calculations a kpts file has to be provided'
		end if
		!
		return
	end subroutine


	subroutine k_space_allocator( H_tb, r_tb, en_k, V_ka, A_ka, Om_kab)
		real(dp),		allocatable			::	en_k(:)
		complex(dp),	allocatable			::	H_tb(:,:,:), r_tb(:,:,:,:), 			&
												A_ka(:,:,:), Om_kab(:,:,:,:),			&
												V_ka(:,:,:)
		!
		allocate(	en_k(						size(H_tb,2)	)	)
		allocate(	V_ka(	3,	size(H_tb,1),	size(H_tb,2)	)	)
		if(	allocated(r_tb)	)	then	
			allocate(	A_ka(	3,		size(r_tb,2),	size(r_tb,3)	)	)
			allocate(	Om_kab(	3,3,	size(r_tb,2),	size(r_tb,3)	)	)
		end if
		!
		return
	end subroutine

end module