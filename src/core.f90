module core
	!this module uses a semiclassic approach to calculate the first ordrer correction
	!	to the polariztion induced by a perturbive magnetic field
	! 	see Niu PRL 112, 166601 (2014)
	!use omp_lib
	use mpi
	use matrix_math,	only:	my_Levi_Civita
	use constants,		only:	dp, aUtoAngstrm, auToTesla,						&			
								mpi_root_id, mpi_id, mpi_nProcs, ierr			
	use statistics,		only:	fd_get_N_el
	use input_paras,	only:	a_latt, 										&
								valence_bands, 									&
								seed_name,										&
								kubo_tol,										&
								hw, eFermi, T_kelvin, i_eta_smr, 				&
								unit_vol
	use file_io,		only:	mpi_read_tb_basis,								&
								write_mep_tensors,								&
								write_ahc_tensor,								&
								write_opt_tensors
	use k_space,		only:	get_recip_latt, get_mp_grid, 					&
								get_rel_kpt,									&
								normalize_k_int
	use wann_interp,	only:	get_wann_interp
	!
	use mep_niu,		only:	mep_niu_CS,	mep_niu_IC, mep_niu_LC
	use kubo,			only:	kubo_ahc_tens,									&
								kubo_opt_tens
		
	implicit none



	private
	public ::			core_worker
	!
	save
	integer									::		num_kpts
contains



!TODO CHECK INDEXING OF VELO ACONN FCURV AFTER THEY HAVE BEEN CALCULATED


!public
	subroutine	core_worker()
		!	interpolate the linear magnetoelectric response tensor alpha
		!
		!	lc	: local current contribution
		!	ic	: itinerant current contribution
		!	cs	: chern - simons term	
		!
		real(dp)							::	N_el_k, N_el_loc, N_el_glob,			&
												kpt(3),	recip_latt(3,3),				&
												!local sum targets:
												mep_tens_ic_loc(	3,3),	 			&
												mep_tens_lc_loc(	3,3),				&
												mep_tens_cs_loc(	3,3),				&
												kubo_ahc_loc(		3,3),				&
												!mpi sum targets:
												mep_tens_ic_glob(	3,3),				&
												mep_tens_lc_glob(	3,3),				&
												mep_tens_cs_glob(	3,3),				&
												kubo_ahc_glob(		3,3)						
												!
		complex(dp)							::	kubo_opt_s_loc(		3,3),				&
												kubo_opt_a_loc(		3,3),				&												
												kubo_opt_s_glob(	3,3),				&
												kubo_opt_a_glob(	3,3),				&
												tempS(3,3), tempA(3,3)
												!
		integer								::	kix, kiy, kiz, ki, 						&
												mp_grid(3), n_ki_loc, n_ki_glob
		complex(dp),	allocatable			::	H_tb(:,:,:), r_tb(:,:,:,:), 			&
												A_ka(:,:,:), Om_kab(:,:,:,:),			&
												V_ka(:,:,:)
		real(dp),		allocatable			::	en_k(:), R_vect(:,:)
		!
		!
		!	k-space
		mp_grid				=	get_mp_grid()
		n_ki_loc			= 	0
		num_kpts			= 	mp_grid(1)*mp_grid(2)*mp_grid(3)
		recip_latt			=	get_recip_latt()
		!
		!read real & allocate k-space
		call mpi_read_tb_basis(	seed_name, R_vect,		 	H_tb, r_tb					)
		call kspace_allocator(		H_tb, r_tb, 			en_k, V_ka, A_ka, Om_kab	)
		!
		!	electron count
		N_el_k		=	0.0_dp
		N_el_loc	=	0.0_dp
		N_el_glob	=	0.0_dp
		!
		!	MagnetoElectric Polarizability
		mep_tens_ic_loc		=	0.0_dp
		mep_tens_lc_loc		=	0.0_dp
		mep_tens_cs_loc		=	0.0_dp
		mep_tens_ic_glob	=	0.0_dp
		mep_tens_lc_glob	=	0.0_dp
		mep_tens_cs_glob	=	0.0_dp
		! 	Anomalous Hall Conducitivity
		kubo_ahc_loc		=	0.0_dp
		kubo_ahc_glob		=	0.0_dp
		!	Optical Conductivity
		kubo_opt_a_loc		=	0.0_dp
		kubo_opt_s_loc		=	0.0_dp
		kubo_opt_a_glob		=	0.0_dp
		kubo_opt_s_glob		=	0.0_dp
		!
		!	
		write(*,'(a,i3,a,i4,a)')		"[#",mpi_id,"; berry_worker]: I start interpolating now (nValence=",valence_bands,")."
		do kiz = 1, mp_grid(3)
			do kiy = 1, mp_grid(2)
				do kix = 1, mp_grid(1)
					!
					!
					ki	=	get_rel_kpt(kix, kiy, kiz, kpt)
					if( mpi_ki_selector(ki, num_kpts)	) then
						!
						!
						call get_wann_interp(H_tb, r_tb, a_latt, recip_latt, R_vect, kpt(:), 	en_k, V_ka, A_ka, Om_kab )
						!
						!MEP
						mep_tens_ic_loc	=	mep_tens_ic_loc + 	mep_niu_IC(V_ka, en_k)		!	itinerant		(Kubo)
						mep_tens_lc_loc	=	mep_tens_lc_loc + 	mep_niu_LC(V_ka, en_k)		!	local			(Kubo)
						mep_tens_cs_loc =	mep_tens_cs_loc + 	mep_niu_CS(A_ka, Om_kab)		!	chern simons	(geometrical)
						!
						!AHC
						kubo_ahc_loc	=	kubo_ahc_loc	+ 	kubo_ahc_tens(en_k,	Om_kab, eFermi, T_kelvin, unit_vol)
						!
						!OPT
						!hw, vol, e_fermi, T_kelvin, i_eta_smr , en_k, A_ka, opt_symm, opt_asym
						call kubo_opt_tens(hw, unit_vol, eFermi, T_kelvin, i_eta_smr, en_k, A_ka, 	tempS, tempA)
						kubo_opt_s_loc	=	kubo_opt_s_loc	+	tempS							
						kubo_opt_a_loc	=	kubo_opt_a_loc	+	tempA
						!
						!COUNT ELECTRONS
						N_el_k	=	get_N_el(en_k, eFermi, T_kelvin)

						!
						n_ki_loc = n_ki_loc + 1
					end if
					!
					!
				end do
			end do
		end do
		write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; berry_worker]: finished interpolating ",n_ki_loc," kpts"
		!
		!
		!sum mep over global kpts
		call MPI_REDUCE(	n_ki_loc,				n_ki_glob,			1,		MPI_INTEGER,		MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	mep_tens_ic_loc,	mep_tens_ic_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		, 	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	mep_tens_lc_loc,	mep_tens_lc_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		, 	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	mep_tens_cs_loc,	mep_tens_cs_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	kubo_ahc_loc,		kubo_ahc_glob,			9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	kubo_opt_s_loc,		kubo_opt_s_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD,	ierr)
		call MPI_REDUCE(	kubo_opt_a_loc,		kubo_opt_a_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD,	ierr)

		!
		!write some files
		if(mpi_id == mpi_root_id) 	then
			write(*,*)	"--------------------------------------------------------------------------------------------------------"	
			if(	n_ki_glob	/=	num_kpts	) stop "[berry_worker]:	ERROR n_ki_glob is not equal to the given mp_grid"
			!
			call normalize_k_int(mep_tens_ic_glob)
			call normalize_k_int(mep_tens_lc_glob)
			call normalize_k_int(mep_tens_cs_glob)
			!
			kubo_ahc_glob	=	kubo_ahc_glob		/	real(n_ki_glob,dp)
			kubo_opt_s_glob	=	kubo_opt_s_glob		/	real(n_ki_glob,dp)			
			kubo_opt_a_glob	=	kubo_opt_a_glob		/	real(n_ki_glob,dp)

			!
			call write_mep_tensors(mep_tens_ic_glob, mep_tens_lc_glob, mep_tens_cs_glob)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; berry_worker]: calculated MEP tensor on ",n_ki_glob," kpts"
			call write_ahc_tensor(kubo_ahc_glob)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; berry_worker]: calculated AHC tensor on ",n_ki_glob," kpts"
			call write_opt_tensors(kubo_opt_s_glob, kubo_opt_a_glob)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; berry_worker]: calculated OPT tensor on ",n_ki_glob," kpts"
		end if
		!
		!
		return
	end subroutine









!private:
	!
	!
!HELPERS
	subroutine kspace_allocator(H_tb, r_tb, en_k, V_ka, A_ka, Om_kab)
		complex(dp),	allocatable, intent(inout)	::		H_tb(:,:,:), r_tb(:,:,:,:),					& 
															V_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp),		allocatable, intent(inout)	::		en_k(:)
		!
		!allocate k-space
		allocate(	en_k(						size(H_tb,2)	)	)
		allocate(	V_ka(	3,	size(H_tb,1),	size(H_tb,2)	)	)
		!
		if(	allocated(r_tb)	)	then	
			allocate(	A_ka(	  3,		size(r_tb,2),	size(r_tb,3)	)	)
			allocate(	Om_kab(	3,	3,		size(r_tb,2),	size(r_tb,3)	)	)
			write(*,'(a,i3,a)')		"[#",mpi_id,"; berry_worker]: will use position operator"	
		end if
		!
		return
	end subroutine


	logical pure function mpi_ki_selector(ki_request, num_kpts)
		integer,		intent(in)		::		ki_request, num_kpts 
		integer							::		ki_todo
		!
		mpi_ki_selector = .false.
		!
		loop_todos: do ki_todo = mpi_id +1, num_kpts, mpi_nProcs
			if(	ki_request == ki_todo) then
				mpi_ki_selector = .true.
				exit loop_todos
			end if
		end do loop_todos
		!
		return
	end function






end module core


