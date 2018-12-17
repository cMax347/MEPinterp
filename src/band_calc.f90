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
	use wann_interp,	only:			get_wann_interp, W_to_H_gaugeTRAFO									

	implicit none


	private
	public		::						band_worker


contains

	subroutine band_worker()
		real(dp),		allocatable			::	rel_kpts(:,:), en_k(:)
		real(dp)							::	recip_latt(3,3)
		integer								::	num_kpts, num_wann, ki, k_per_mpi
		complex(dp),	allocatable			::	H_k(:,:),	&
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
				call get_ham(rel_kpts(:,ki),	H_k,	V_ka	)
				!
				!
				if(do_gauge_trafo)	call do_jans_gauge_trafo(ki, H_k, V_ka, en_k)
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



	subroutine do_jans_gauge_trafo(kpt_idx, H_k, V_ka, en_k)
		integer,			intent(in)		::	kpt_idx
		complex(dp),		intent(in)		::	H_k(:,:)
		complex(dp),		intent(inout)	::	V_ka(:,:,:)
		real(dp),			intent(out)		::	en_k(:)
		integer, 			allocatable 	::	iwork(:),ifail(:)
      	real(dp),			allocatable 	::	rwork(:)
		complex(dp),		allocatable		::	work(:), z(:,:), A_ka(:,:,:), Om_kab(:,:,:,:)
		integer								::	num_wann, lwork,lrwork,liwork,info,nev, i, x
		real(dp) 							::	vl,vu,abstol
		!
		!	JANS SOURCE CODE:
		!
		!       call init_ham( kpts(1:3,k),num_wann,ham,
		!     >                vW(:,:,1),vW(:,:,2),vW(:,:,3) )
		!
		!       call zheevx('V','A','U',num_wann,ham,num_wann,
		!     >             vl,vu,1, num_wann,abstol,nev,
		!     >             eig(:,k),z(:,:),num_wann,work,lwork,
		!     >             rwork,iwork,ifail,info)
		!       if(info.ne.0) stop 'zheevx'
		!
		!       do i=1,3
		!        vH(:,:,i) = matmul(matmul(conjg(transpose(z(:,:))),
		!     >                   vW(:,:,i)),z(:,:))
		!       enddo
		!
		abstol 		= 	0.0 !2.0*tiny(abstol)
		num_wann	=	8
      	lwork		=	12*num_wann
     	lrwork		=	17*num_wann
      	liwork		=	15*num_wann
      	!
      	allocate( z(num_wann, num_wann))
		allocate( rwork(lrwork),work(lwork),iwork(liwork) )
		allocate( ifail(num_wann) )
		!
		!
		!!!	GET EIGENVECTORS
		!call zheevx('V','A','U',num_wann,H_K,num_wann,				&
		!	             vl,vu,1, num_wann,abstol,nev,				&
		!	             en_k(:),z(:,:),num_wann,work,lwork,		&
		!	             rwork,iwork,ifail,info						&
		!	        )
		!if(info /= 0 ) stop 'zheevx'
		
		!z(:,:)	=	H_k(:,:)
		!call zheevd_wrapper(z, en_k)


		!W_to_H_gaugeTRAFO(e_k, U_k, H_ka, A_ka, Om_kab)
		z(:,:)	=	H_k(:,:)
		call zheevd_wrapper(z, en_k)



		!call W_to_H_gaugeTRAFO(en_k, z(:,:), V_ka(:,:,:), A_ka, Om_kab)




		!	PERFORM GAUGE TRAFO
		do i = 1, 3
			V_ka(i,:,:)		=	matmul(		matmul( conjg(transpose(z(:,:))), V_ka(i,:,:)), 		z(:,:)		)
		end do
		write(*,'(a,i7)')	"[do_jans_gauge_trafo]: finished gauging velos at #kpt",kpt_idx	
		!
		return
	end subroutine


end module