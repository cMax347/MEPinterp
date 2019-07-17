module input_paras
	use m_config
#ifdef USE_MPI
	use mpi
#endif
	use matrix_math,				only:		crossP
	use constants,					only:		dp, fp_acc, pi_dp, aUtoEv, kBoltz_Eh_K
	use mpi_community,				only:		mpi_id, mpi_root_id, mpi_nProcs, ierr
	use k_space,					only:		set_recip_latt, 	&
												set_mp_grid,		& 
												set_kspace_cutoff, 	&
												print_kSpace_info




	implicit none

	private
	public								::		&
												!routines		
												init_parameters, 	my_mkdir,						 		 		&
												!dirs
												w90_dir, 															&
												raw_dir,															&
												out_dir,															&
												velo_out_dir,														&
												!jobs
												plot_bands,															&
												use_mpi,															&
												debug_mode,															&
												use_kspace_ham, kspace_ham_id,										&
												use_cart_velo,														&
												do_gauge_trafo,														&
												do_write_velo,														&
												use_R_float,														&
												do_write_mep_bands,													&
												do_mep, do_ahc, do_kubo, do_opt, do_photoC, do_bcd_photo, do_gyro,	&
												!atoms
												wf_centers,															&
												!vars
												seed_name,	valence_bands,											&
												a_latt, kubo_tol, unit_vol,											&
												k_cutoff,															&
												N_hw, hw_min, hw_max,												&
												N_eF, eF_min, eF_max, T_kelvin,										&
												N_eta_smr ,eta_smr_min, eta_smr_max



	
	save

	
	integer						::	valence_bands
	character(len=:), allocatable::	seed_name
	character(len=9)			::	w90_dir	="w90files/"
	character(len=4)			::	raw_dir ="raw/"
	character(len=4)			::	out_dir	="out/"
	character(len=10)			:: 	velo_out_dir
	logical						::	plot_bands, 					&
									use_kspace_ham,					&
									use_cart_velo,					&
									do_gauge_trafo, 				&
									use_R_float,					&
									do_write_velo,					&
									do_write_mep_bands,				&
									debug_mode,	use_mpi,			&
									do_mep, do_ahc, do_kubo, 		&
									do_opt, do_photoC, do_bcd_photo,& 
									do_gyro
	integer						::	N_wf, N_eta_smr, N_eF, N_hw, kspace_ham_id
	real(dp)					::	a_latt(3,3), a0, unit_vol,		&
									kubo_tol, hw_min, hw_max,		&
									k_cutoff=1.0_dp,				&
									eF_min, eF_max, 				&
									eta_smr_min, eta_smr_max,		&
									T_kelvin			
	real(dp),	allocatable		::	wf_centers(:,:)




	contains




!public
	logical function init_parameters()
		type(CFG_t) 			:: 	my_cfg
		real(dp)				::	a1(3), a2(3), a3(3)
		integer					::	mp_grid(3)
		logical					::	input_exist
		character(len=132)		:: 	long_seed_name
		!
		use_mpi	= .false.
#ifdef USE_MPI
		use_mpi = .true.
#endif
		init_parameters	=	.false.
		input_exist		=	.false.
		!
		if( 	.not. use_mpi 		.and.		 mpi_id	/= 0			)	then
			 write(*,'(a,i7.7,a)')		'[#',mpi_id,';init_parameters]:	hello, I am an unexpected MPI thread !!!1!1!!!1!!1 '
		end if
		!
		velo_out_dir	=	out_dir//"/velo/"
		kspace_ham_id	=	-1
		!ROOT READ
		if(mpi_id == mpi_root_id) then
			inquire(file="./input.cfg",exist=input_exist)
			!
			if( input_exist)	then		
				!OPEN FILE
				call CFG_read_file(my_cfg,"./input.cfg")
				!
				![methods]
				call CFG_add_get(my_cfg,	"jobs%plot_bands"				,	plot_bands			,	"if true do a bandstructure run"		)
				call CFG_add_get(my_cfg,	"jobs%debug_mode"				,	debug_mode			,	"switch aditional debug tests in code"	)
				call CFG_add_get(my_cfg,	"jobs%R_vect_float"				,	use_R_float			,	"the R_cell vector is now real (else: integer)")
				call CFG_add_get(my_cfg,	"jobs%do_write_velo"			,	do_write_velo		,	"write formatted velocity files at each kpt")
				call CFG_add_get(my_cfg,	"jobs%do_mep"					,	do_mep				,	"switch (on/off) this response tens calc")
				call CFG_add_get(my_cfg,	"jobs%do_kubo"					,	do_kubo				,	"switch (on/off) this response tens calc")
				call CFG_add_get(my_cfg,	"jobs%do_ahc"					,	do_ahc				,	"switch (on/off) this response tens calc")
				call CFG_add_get(my_cfg,	"jobs%do_opt"					,	do_opt				,	"switch (on/off) this response tens calc")
				call CFG_add_get(my_cfg,	"jobs%do_photoC"				,	do_photoC			,	"switch (on/off) this response tens calc")				
				call CFG_add_get(my_cfg,	"jobs%do_bcd_photo"				,	do_bcd_photo		,	"switch (on/off) this response tens calc")
				call CFG_add_get(my_cfg,	"jobs%do_gyro"					,	do_gyro				,	"switch (on/off) this response tens calc")
				!~~~~~~~~~~~~
				!
				![unitCell]
				call CFG_add_get(my_cfg,	"unitCell%a1"      				,	a1(1:3)  	  		,	"a_x lattice vector	(Bohr)"				)
				call CFG_add_get(my_cfg,	"unitCell%a2"      				,	a2(1:3)  	  		,	"a_y lattice vector	(Bohr)"				)
				call CFG_add_get(my_cfg,	"unitCell%a3"      				,	a3(1:3)  	  		,	"a_z lattice vector	(Bohr)"				)
				call CFG_add_get(my_cfg,	"unitCell%a0"					,	a0					,	"lattice scaling factor "				)
				!~~~~~~~~~~~~
				!
				![wannBase]
				call CFG_add_get(my_cfg,	"wannBase%seed_name"			,	long_seed_name		,	"seed name of the TB files"				)
				call CFG_add_get(my_cfg,	"wannBase%N_wf"					,	N_wf				,	"number of WFs specified in input"		)
				if(	N_wf > 0) then 
					allocate(	wf_centers(		3,	N_wf)	)
					call CFG_add_get(my_cfg,	"wannBase%wf_centers_x"		,	wf_centers(1,:)			,	"array of x coord of relative pos"	)
					call CFG_add_get(my_cfg,	"wannBase%wf_centers_y"		,	wf_centers(2,:)			,	"array of y coord of relative pos"	)
					call CFG_add_get(my_cfg,	"wannBase%wf_centers_z"		,	wf_centers(3,:)			,	"array of z coord of relative pos"	)
				end if
				call CFG_add_get(my_cfg,	"wannBase%use_kspace_ham"		,	use_kspace_ham		,	"swith for using model ham setup in k-space")
 				if(use_kspace_ham)	&
 					call CFG_add_get(my_cfg,	"wannBase%kspace_ham_id"	,	kspace_ham_id		,	"IDs the model to use (0:rashba)	"	)
				call CFG_add_get(my_cfg,	"wannBase%k_cutoff"				,	k_cutoff			,	"cutoff for internal k values"			)
				!~~~~~~~~~~~~
				!
				![wannInterp]
				call CFG_add_get(my_cfg,	"wannInterp%use_cart_velo"		,	use_cart_velo		,	"use cartesian instead of internal units")
				call CFG_add_get(my_cfg,	"wannInterp%doGaugeTrafo"		,	do_gauge_trafo		,	"switch (W)->(H) gauge trafo"			)
				call CFG_add_get(my_cfg,	"wannInterp%mp_grid"			,	mp_grid(1:3)		,	"interpolation k-mesh"					)
				!~~~~~~~~~~~~
				!
				![Fermi]
				call CFG_add_get(my_cfg,	"Fermi%N_eF"					,	N_eF				,	"number of fermi energys to test"		)
				call CFG_add_get(my_cfg,	"Fermi%eF_min"					,	eF_min				,	"minimum fermi energy( in eV)"			)
				call CFG_add_get(my_cfg,	"Fermi%eF_max"					,	eF_max				,	"maximum fermi energy( in eV)"			)
				call CFG_add_get(my_cfg,	"Fermi%Tkelvin"					,	T_kelvin			,	"Temperature"							)				
				call CFG_add_get(my_cfg,	"Fermi%N_eta_smr"				,	N_eta_smr			,	"number of electr smr to be probed"		)
				call CFG_add_get(my_cfg,	"Fermi%eta_smr_min"				,	eta_smr_min			,	"min smearing (eV) applied"					)
				call CFG_add_get(my_cfg,	"Fermi%eta_smr_max"				,	eta_smr_max			,	"max smearing (eV) applied"					)
				call CFG_add_get(my_cfg,	"Fermi%kuboTol"					,	kubo_tol			,	"numerical tolearnce for KUBO formulas"	)
				!~~~~~~~~~~~~
				!
				![mep]
				call CFG_add_get(my_cfg,	"MEP%valence_bands"				,	valence_bands		,	"number of valence_bands"				)
				call CFG_add_get(my_cfg,	"MEP%do_write_mep_bands"		,	do_write_mep_bands	,	"write mep tensor band resolved"		)
				!~~~~~~~~~~~~
				
				!
				![Laser]
				call CFG_add_get(my_cfg,	"Laser%N_hw"					,	N_hw					,	"points to probe in interval"		)
				call CFG_add_get(my_cfg,	"Laser%hw_min"					,	hw_min					,	"min energy of incoming light"		)
				call CFG_add_get(my_cfg,	"Laser%hw_max"					,	hw_max					,	"max energy of incoming light"		)
				!~~~~~~~~~~~~
				!~~~~~~~~~~~~
				!~~~~~~~~~~~~
				!
				!	lattice setup
				a_latt(1,1:3)	= a1(1:3)
				a_latt(2,1:3)	= a2(1:3)
				a_latt(3,1:3)	= a3(1:3)
				a_latt			=	a0 * a_latt
				!
				! 	unit conversion
				hw_min			=	hw_min 		/ 	aUtoEv
				hw_max			=	hw_max 		/ 	aUtoEv
				!
				N_hw			=	max(1,N_hw)
				N_ef			=	max(1,N_ef)
				N_eta_smr		=	max(1,N_eta_smr)
				!
				eF_min			= 	eF_min		/	aUtoEv
				eF_max			=	eF_max		/	aUtoEv
				!
				eta_smr_min		=	eta_smr_min	/	aUtoEv
				eta_smr_max		=	eta_smr_max	/	aUtoEv


				if(.not.use_kspace_ham)then	
					if(abs(k_cutoff-1.0_dp)>1e-3_dp) then
						write(*,*)	"[#",mpi_id,";init_parameters]: WARNING k_cutoff set to 1.0 (only allowed if use_kspace_ham is true)"
					end if
					k_cutoff = 1.0_dp
				end if
				!
				!
				!	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
				write(*,*)					"*"
				write(*,*)					"*"
				write(*,*)					"parallelization with #",mpi_nProcs," MPI threads"
				write(*,*)					""
				write(*,*)					"**********************init_parameter interpretation******************************"
				write(*,'(a,i7.7,a)')			"[#",mpi_id,";init_parameters]: input interpretation:"
				!call CFG_write(my_cfg,	"stdout", hide_unused=.false.)	!writes all input vars to stdout
				write(*,*)					"~~"
				if(	T_kelvin < 1e-2_dp)	 then
					write(*,'(a)')			"	WARNING , too small temperature value (<1e-2). WARNING Fermi Dirac will assume T=0 (stepfunction)"
				else
					write(*,*)				"thermal smearing :"	,	kBoltz_Eh_K *	T_kelvin / aUtoEv		,	" (eV)"
				end if

				!	------------------------------------------
				!	------------------------------------------					
				!	------------------------------------------
				
				
				write(*,*)					"*********************************************************************************"		
				!call CFG_write(my_cfg, './input.log', hide_unused=.true.)
				write(*,'(a,i7.7,a)')			"[#",mpi_id,";init_parameters]: wrote input log file to input.log"
				write(*,*)		""
				!
				!make the output folder
				write(*,*)	"*"
				write(*,*)	"----------------------MAKE TARGET DIRECTORIES-------------------------"
				write(*,'(a,i7.7,a)')			"[#",mpi_id,";init_parameters]: start target mkdir..."
				call my_mkdir(out_dir)
				call my_mkdir(raw_dir)
				if(		 	do_write_velo		)	call my_mkdir(velo_out_dir)	
				write(*,'(a,i7.7,a)')			"[#",mpi_id,";init_parameters]: ...all required directories created"
			else
				write(*,'(a,i7.7,a)')			"[#",mpi_id,";init_parameters]: could not find input file"
				stop "please provide a input.cfg file"
			end if
			write(*,*)	"*"
			write(*,*)	"----------------------K-SPACE SETUP-------------------------"
			write(*,'(a,i7.7,a)')			"[#",mpi_id,";init_parameters]: now bcast the input parameters and setup k-space ..."
			write(*,*)	""
		end if
		!
		!
		if(use_mpi) then
			call MPI_BCAST(		input_exist		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
		endif
		!
		!
		if( input_exist) then
			if(use_mpi) then
				!ROOT BCAST
				![FLAGS]
				call MPI_BCAST(		plot_bands		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		debug_mode		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		use_kspace_ham	,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		use_cart_velo	,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		do_gauge_trafo	,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		do_write_velo	,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_write_mep_bands,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		use_R_float		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_mep 			,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_kubo 		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_ahc 			,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_opt 			,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_photoC		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_bcd_photo	,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		do_gyro 		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				![SYSTEM]
				call MPI_BCAST(		kspace_ham_id	,			1			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		a_latt			,			9			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		valence_bands	,			1			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		long_seed_name	,	len(long_seed_name)	,		MPI_CHARACTER		,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		mp_grid			,			3			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		k_cutoff		,			1			,		MPI_DOUBLE_PRECISION,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				![ATOMS]
				call MPI_BCAST(		N_wf			,			1			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				if(	N_wf > 0) then
					if(	mpi_id /= mpi_root_id )		allocate(	wf_centers(3,N_wf)	)
					call MPI_BCAST(		wf_centers			,		3* N_wf			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				end if
				![KUBO]
				call MPI_BCAST(		hw_min			,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)			
				call MPI_BCAST(		hw_max			,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		n_hw			,			1			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				![FERMI]
				call MPI_BCAST(		kubo_tol		,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD,	ierr)	
				call MPI_BCAST(		N_eF			,			1			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		eF_min			,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		eF_max			,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		T_kelvin		,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		N_eta_smr		,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		eta_smr_min		,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		eta_smr_max		,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
			end if
			!
			!TRIM SEEDNAME
			allocate(	character(len=len(trim(long_seed_name)))	::	seed_name		)
			seed_name	=	trim(long_seed_name)
			!
			!UNIT CELL VOLUME
			a1			=	a_latt(1,:)
			a2			=	a_latt(2,:)
			a3			=	a_latt(3,:)
			unit_vol	=	dot_product(	crossP(a1, a2)	,	a3		)
			!
			!SETUP K-SPACE
			call set_recip_latt(a_latt)
			call set_mp_grid(mp_grid)
			call set_kspace_cutoff(k_cutoff)
			call print_kSpace_info()
		end if
		!
		init_parameters	=	input_exist
		!
		!
		return 
	end function





	subroutine my_mkdir(dir)
		character(len=*)			::	dir
		!logical						::	dir_exists
		character(len=11)		::	mkdir="mkdir -p ./"	!to use with system(mkdir//$dir_path) 	
		!
		!inquire(directory=dir, exist=dir_exists)
		!if( .not. dir_exists )	then
			call system(mkdir//dir)
			write(*,'(a,i7.7,a,a)')	"[#",mpi_id,"; init_parameters]: (fake) created directory ",dir
		!end if
		!
		return
	end subroutine




end module input_paras