module band_calc
!
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
#ifdef USE_MPI
	use mpi
#endif
	use constants,		only:			dp
	use mpi_comm,		only:			mpi_root_id, mpi_id, mpi_nProcs, ierr

	use input_paras,	only:			use_mpi, seed_name, a_latt, 					&
										do_gauge_trafo, do_write_velo
	use k_space,		only:			get_recip_latt										
	use wrapper_3q,		only:			get_ham
	use matrix_math,	only:			zheevd_wrapper, zheevr_wrapper
	use file_io,		only:			read_kptsgen_pl_file,							&
										write_en_binary, 								&
										write_en_global,								&
										write_velo
	use wann_interp,	only:			get_wann_interp									

	implicit none


	private
	public		::						band_worker


contains

	subroutine band_worker()
		real(dp),		allocatable			::	rel_kpts(:,:), en_k(:)
		real(dp)							::	recip_latt(3,3)
		integer								::	num_kpts, ki, k_per_mpi
		complex(dp),	allocatable			::	H_tb(:,:,:), r_tb(:,:,:,:), 				&
												A_ka(:,:,:), Om_kab(:,:,:,:),				&
												V_ka(:,:,:)					
		!
		if(mpi_id==mpi_root_id)	then
			write(*,*)	"----------------------------------------------------------------"
			write(*,*)	"----------------------------------------------------------------"
			write(*,*)	"----------------------------------------------------------------"
			write(*,*)	"*"
			write(*,*)	"*"
			write(*,*)	"***^^^^	-	BANDSTRUCTURE MODE	-	^^^^***"
			write(*,*)	"*"
			write(*,*)	"*"
			write(*,*)	"----------------------------------------------------------------"
			write(*,*)	"----------------------------------------------------------------"
			write(*,*)	"----------------------------------------------------------------"
		end if
		!
		if( read_kptsgen_pl_file(rel_kpts)	) then
			!
			!	get k-space
			recip_latt	= get_recip_latt()
			num_kpts	= size(rel_kpts,2)
			k_per_mpi	= 0
			!
			!	allocate
			num_wann	=	8
			allocate(	en_k(		num_wann			)	)
			allocate(	H_k(	num_wann,	num_wann	)	)
			!
			!
			!	do the work
			if(mpi_nProcs>1)		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
			write(*,'(a,i3,a,a,a)')	'[#',mpi_id,';band_worker/',cTIME(time()),	']:	start interpolating...'
			!
			!
			do ki = mpi_id + 1, num_kpts,	mpi_nProcs
				!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
				!			ONLY GET HAM															 |
				!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				call get_ham(rel_kpts(:,ki),	H_k,	v_dummy	)
				call zheevd_wrapper(H_k, en_k)
				!
				call write_en_binary(ki,en_k)
				if(do_write_velo)	call write_velo(ki,V_ka)

				k_per_mpi	= k_per_mpi + 1
			end do
			write(*,'(a,i3,a,a,a,i10,a,i10,a)')	'[#',mpi_id,';band_worker/',cTIME(time()),				&
													']:	...finished, interpolated ',k_per_mpi,' of ',num_kpts,' kpts'
			!
			!
			!	write the results
			if(mpi_nProcs>1)	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
			if(mpi_id == mpi_root_id)	then
				call write_en_global(num_wann,	rel_kpts)
				write(*,*)'---------------------------------------------------------------------------------------------'
			end if
		else
			stop 'for bandstructure calculations a kpts file has to be provided'
		end if
		!
		return
	end subroutine


end module