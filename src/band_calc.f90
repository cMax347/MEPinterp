module band_calc
	use mpi
	use constants,		only:			dp, mpi_root_id, mpi_id, mpi_nProcs, ierr
	use input_paras,	only:			seed_name, a_latt
	use k_space,		only:			get_recip_latt										
	use file_io,		only:			read_kptsgen_pl_file,							&
										mpi_read_tb_basis,								&
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
		!
		if( read_kptsgen_pl_file(rel_kpts)	) then
			recip_latt	= get_recip_latt()
			num_kpts	= size(rel_kpts,2)
			k_per_mpi	= 0
			!
			!	get the data
			call mpi_read_tb_basis(seed_name, R_vect, H_tb, r_tb)
			!
			allocate(	en_k(						size(H_tb,2)	)	)
			allocate(	V_ka(	3,	size(H_tb,1),	size(H_tb,2)	)	)
			if(	allocated(r_tb)	)	then	
				allocate(	A_ka(	3,		size(r_tb,2),	size(r_tb,3)	)	)
				allocate(	Om_kab(	3,3,	size(r_tb,2),	size(r_tb,3)	)	)
				!write(*,'(a,i3,a)')		"[#",mpi_id,"; mep_interp]: will use position operator"
			end if
			!
			!	do the work
			do ki = mpi_id + 1, num_kpts,	mpi_nProcs
				call get_wann_interp(H_tb, r_tb, a_latt, recip_latt, R_vect, rel_kpts(:,ki), 	en_k, V_ka, A_ka, Om_kab )
				call write_en_binary(ki,en_k)
				k_per_mpi	= k_per_mpi + 1
			end do
			write(*,'(a,i10,a,i10,a)')	'[band_worker]:	interpolated ',k_per_mpi,' of ',num_kpts,' kpts'
			!
			!
			!	write the results
			call MPI_BARRIER(MPI_COMM_WORLD, ierr)
			write(*,*) 'beyond the barrier'
			if(mpi_id == mpi_root_id)	call write_en_global(rel_kpts)
		else
			stop 'for bandstructure calculations a kpts file has to be provided'
		end if
		!
		return
	end subroutine




end module