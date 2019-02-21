module band_calc
!
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
#ifdef USE_MPI
	use mpi
#endif
	use constants,		only:			dp
	use mpi_community,	only:			mpi_root_id, mpi_id, mpi_nProcs, ierr
	use input_paras,	only:			use_mpi, seed_name, a_latt, 					&
										do_gauge_trafo, 								& 
										atPos,											&
										do_write_velo
	use file_io,		only:			read_kptsgen_pl_file,							&
										read_tb_basis,									&
										write_en_binary, 								&
										write_en_global,								&
										write_velo
	use wann_interp,	only:			get_wann_interp									
	implicit none


	private
	public		::						band_worker


contains

	subroutine band_worker()
		real(dp),		allocatable			::	rel_kpts(:,:), en_k(:), R_vect(:,:)
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
			num_kpts	= size(rel_kpts,2)
			k_per_mpi	= 0
			!
			!	get the data
			call read_tb_basis(seed_name, R_vect, H_tb, r_tb)
			call band_allocator(H_tb, r_tb, en_k, V_ka, A_ka, Om_kab)
			!
			!
			!	do the work
			if(use_mpi)		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
			write(*,'(a,i3,a,a,a)')	'[#',mpi_id,';band_worker/',cTIME(time()),	']:	start interpolating...'
			!
			!
			do ki = mpi_id + 1, num_kpts,	mpi_nProcs
				call get_wann_interp(	do_gauge_trafo, H_tb, r_tb, 								&
										R_vect, 	atPos	,  ki, rel_kpts(:,ki),							& 
										en_k, V_ka, A_ka, Om_kab									&
									)
				call write_en_binary(ki,en_k)
				if(do_write_velo)	call write_velo(ki,V_ka)

				k_per_mpi	= k_per_mpi + 1
			end do
			write(*,'(a,i3,a,a,a,i10,a,i10,a)')	'[#',mpi_id,';band_worker/',cTIME(time()),				&
													']:	...finished, interpolated ',k_per_mpi,' of ',num_kpts,' kpts'
			!
			!
			!	write the results
			if(use_mpi)	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
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


	subroutine band_allocator( H_tb, r_tb, en_k, V_ka, A_ka, Om_kab)
		real(dp),		allocatable, intent(inout)			::	en_k(:)
		complex(dp),	allocatable, intent(inout)			::	H_tb(:,:,:), r_tb(:,:,:,:), 	&
																V_ka(:,:,:),					&
																A_ka(:,:,:), Om_kab(:,:,:,:)									
		!
		allocate(	en_k(size(H_tb,2))		)
		en_k	=	0.0_dp
		!
		if( do_write_velo )		then
			allocate(	V_ka(3,	size(H_tb,1), size(H_tb,2))		)
			V_ka	=	0.0_dp
			write(*,'(a,i3,a,i1,a,i3,a,i3,a)')	"[",mpi_id,"band_allocator]: V_ka allocated with size=	(",&
													size(V_ka,1), "x",size(V_ka,2),"x",size(V_ka,3),")"
		end if
		!
		!
		if(	allocated(r_tb)	)	then	
			allocate(	A_ka(	3,		size(r_tb,2),	size(r_tb,3)	)	)
			allocate(	Om_kab(	3,3,	size(r_tb,2),	size(r_tb,3)	)	)
			A_ka	=	0.0_dp
			Om_kab	=	0.0_dp
		end if
		!
		!
		return
	end subroutine

end module