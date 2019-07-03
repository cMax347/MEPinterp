program test_2Drashba_photoC
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
	use omp_lib
	use mpi
	use m_config
	use constants,				only:		dp
	use mpi_community,			only:		mpi_root_id, mpi_id, mpi_nProcs, ierr
	use input_paras,			only:		init_parameters,	&
											my_mkdir
	use k_space,				only:		get_mp_grid
	use core,					only:		core_worker
	use test_helpers,			only:		my_exit		
	!
	implicit none
	!
	logical							::		test_passed
	!
	
	!
	test_passed	=	.False.
	!MPI INIT
	mpi_root_id = 	0
	mpi_id		=	0
	mpi_nProcs	=	1
	call MPI_INIT( ierr )
    call MPI_COMM_RANK (MPI_COMM_WORLD, 	mpi_id			, ierr)
    call MPI_COMM_SIZE (MPI_COMM_WORLD, 	mpi_nProcs		, ierr)
	!
	!
	if(	mpi_id	==	mpi_root_id	)	then
		write(*,*)					'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
		write(*,*)					"|							2D Rashba TB model test (2nd order photo cond)	|"
		write(*,*)					"|						https://arxiv.org/abs/1811.05803 					|"
		write(*,*)					'-----------------------------------------------------------------------------'
		write(*,*)					'	this test checks:'
		!write(*,*)					'		->	wann interpolation of tight binding Hamiltonian'
		!write(*,*)					'		->	k_space setup & integration'
		!write(*,*)					'		->	MPI parallelization'
		!write(*,*)					'		->	(frequency dependent) AHE'
		!write(*,*)					'		->	tight binding approximation (i.e. wf_centers are used to get FT-phase)'
		!write(*,*)					'		->	curvature interpolation (via zero freq. limit - compare RE{ahcVELO(hw=0)} to ah_tens)'
		write(*,*)					'		-> 	2nd order photo conductivity (Nagaosa formula) in tight binding model'
		write(*,*)					'		-> 	'
		write(*,*)					'	this test DOES NOT CHECK:'
		!write(*,*)					'		 X	position operator interpolation'
		!write(*,*)					'		 X	cartesian velocities in wannier interpolation'
		!write(*,*)					'		 X	any none Hall response tensors'
		write(*,*)					'		 X	influence of position operator on 2nd order photo cond'
		write(*,*)					'		 X	convergence of 2nd order photo cond with respect to number of cond bands included in wann energy window'
		write(*,*)					'		 X	'

		write(*,*)					'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		write(*,*)					'...'
		!				~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		write(*,'(a,i7.7,a,a,a)')	'[#',mpi_id,';test_2Drashba_photoC/',cTIME(time()),']:	hello there'
		call write_input_files()
		write(*,'(a,i7.7,a,a,a)')	'[#',mpi_id,';test_2Drashba_photoC/',cTIME(time()),']:	wrote input files'
	end if
	!
	!
    !READ IN & DO K-SPACE INTEGRATION
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    



    write(*,'(a,i7.7,a,a,a)')	'[#',mpi_id,';test_2Drashba_photoC/',cTIME(time()),']:	start reading input file ...'
    if( init_parameters()	) then
		 write(*,'(a,i7.7,a,a,a)')	'[#',mpi_id,';test_2Drashba_photoC/',cTIME(time()),']:	read input cfg start core worker'
		call core_worker()			
	else
		write(*,'(a,i7.7,a)')		'[#',mpi_id,';test_2Drashba_photoC]: ....input file not found, by'
		call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
	end if
	!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				READ & INTERPRET RESULTS
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		
		!	ToDo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	!FINALIZE
	call MPI_FINALIZE(ierr)
	write(*,*)	'[#',mpi_id,';test_2Drashba_photoC/',cTIME(time()),']:	test was performed on mp_grid ',get_mp_grid()
	write(*,'(a,i7.7,a,a,a)')	'[#',mpi_id,';test_2Drashba_photoC/',cTIME(time()),']:	consider finer mesh for better accuracy'		
	write(*,'(a,i7.7,a,a,a)')	'[#',mpi_id,';test_2Drashba_photoC/',cTIME(time()),']:	all done, by by'	


	if(	test_passed	)	&
		call my_exit(.True.)
	call my_exit(.False.)	






	contains


	subroutine write_input_files()
		! write config file
		integer			::	status_hr, status_npy
		call write_input_cfg_file()
		! write w90 files
			! ->	store w90files in new subfolder of config folder 
			! ->		config/FeMn_tb/			
			!				->	w90files/
			!				->	w90files/
			!				-> ahcVELO.npy
			! ->
			! ->
		call my_mkdir('./w90files/')
		!status_hr	=	system("cp ../../config/FeMn_TB/FeMn_hr.dat	./w90files/")
		!write(*,*)	"[write_input_files]:	copied hr to w90files with status:	",status_hr
		!!
		!status_npy	=	system("cp ../../config/FeMn_TB/ahcVELO.npy	./w90files/ahcVELOsol.npy")
		!write(*,*)	"[write_input_files]:	copied npy file to w90files with status:	",status_npy
		!
		return
	end subroutine



	subroutine write_input_cfg_file()
		!
		real(dp)				::	wf_centers(3,8),			&
									a1(3),a2(3),a3(3),a0
		integer					::	mp_grid(3)
		character(len=13)		::	fname_cfg	= "./input.cfg"
		type(CFG_t)          	:: 	my_cfg
		!
		!
		write(*,*)	"[write_input_cfg_file]:	WARNING ERROR !@!111!! TODO"
		!	LATTICE
		a1(1:3)			=	(/  1.43 , 0.00 ,  0.00	/)
		a2(1:3)			=	(/  0.00 , 1.43 ,  0.00	/)	
		a3(1:3)			=	(/  0.00 , 0.00 ,  1.00	/)
		a0				=	1.0_dp
		!
		!	WF CENTERS
		!		->	spin up centers
		wf_centers(1,1:4)		=	(/	0.0,	0.0,	0.0,	0.0	/)		!	^		
		wf_centers(2,1:4)		=	(/	0.0,	0.0,	0.0,	0.0	/)		!	^		
		wf_centers(3,1:4)		=	(/	0.0,	0.0,	0.0,	0.0	/)		!	^
		!		->	spin dwn centers
		wf_centers(1,1+4:4+4)	=	wf_centers(1,1:4)						!	v						
		wf_centers(2,1+4:4+4)	=	wf_centers(2,1:4)						!	v	
		wf_centers(3,1+4:4+4)	=	wf_centers(3,1:4)						!	v			
		!
		!	NUMERICS
		mp_grid(1:2)	=	10200
		mp_grid(3)		=	1
		!
		!	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		call CFG_add(my_cfg,	"jobs%plot_bands"			,		.False.			,	"if true do a bandstructure run"				)
		call CFG_add(my_cfg,	"jobs%debug_mode"			,		.False.			,	"switch aditional debug tests in code"			)
		call CFG_add(my_cfg,	"jobs%R_vect_float"			,		.True.			,	"the R_cell vector is now real (else: integer)"	)
		call CFG_add(my_cfg,	"jobs%do_write_velo"		,		.False.			,	"write formatted velocity files at each kpt"	)
		call CFG_add(my_cfg,	"jobs%do_mep"				,		.False.			,	"switch (on/off) this response tens calc"		)
		call CFG_add(my_cfg,	"jobs%do_kubo"				,		.False.			,	"switch (on/off) this response tens calc"		)
		call CFG_add(my_cfg,	"jobs%do_ahc"				,		.False.			,	"switch (on/off) this response tens calc"		)
		call CFG_add(my_cfg,	"jobs%do_opt"				,		.False.			,	"switch (on/off) this response tens calc"		)
		call CFG_add(my_cfg,	"jobs%do_photoC"			,		.True.			,	"switch (on/off) this response tens calc"		)
		call CFG_add(my_cfg,	"jobs%do_gyro"				,		.False.			,	"switch (on/off) this response tens calc"		)
		!~~~~~~~~~~~~
		!
		![unitCell]
		call CFG_add(my_cfg,	"unitCell%a1"      			,		a1(1:3)  	  	,	"a_x lattice vector (bohr)"							)
		call CFG_add(my_cfg,	"unitCell%a2"      			,		a2(1:3)  	  	,	"a_y lattice vector (bohr)"							)
		call CFG_add(my_cfg,	"unitCell%a3"      			,		a3(1:3)  	  	,	"a_z lattice vector (bohr)"							)
		call CFG_add(my_cfg,	"unitCell%a0"				,		a0				,	"lattice scaling factor "						)
		!~~~~~~~~~~~~
		!
		![wannBase]
		call CFG_add(my_cfg,	"wannBase%seed_name"		,		"r2d"			,	"seed name of the TB files"						)
		call CFG_add(my_cfg,	"wannBase%N_wf"				,		0				,	"number of WFs specified in input"				)
		call CFG_add(my_cfg,	"wannBase%wf_centers_x"		,	wf_centers(1,:)		,	"array of x coord of relative pos"				)
		call CFG_add(my_cfg,	"wannBase%wf_centers_y"		,	wf_centers(2,:)		,	"array of y coord of relative pos"				)
		call CFG_add(my_cfg,	"wannBase%wf_centers_z"		,	wf_centers(3,:)		,	"array of z coord of relative pos"				)
		call CFG_add(my_cfg,	"wannBase%use_kspace_ham"	,		.True.			,	"swith for using model ham setup in k-space"	)
		call CFG_add(my_cfg,	"wannBase%kspace_ham_id"	,		0				,	"use id=0 for the 2D rasbha model with exchange")
		call CFG_add(my_cfg,	"wannBase%k_cutoff"			,		0.7_dp				,	"k_space cutoff, required to get correct units in open sys")
		!~~~~~~~~~~~~
		!
		![wannInterp]
		call CFG_add(my_cfg,	"wannInterp%use_cart_velo"		,	.False.			,	"use cartesian instead of internal units"		)
		call CFG_add(my_cfg,	"wannInterp%doGaugeTrafo"		,	.True.			,	"switch (W)->(H) gauge trafo"					)
		call CFG_add(my_cfg,	"wannInterp%mp_grid"			,	mp_grid(1:3)	,	"interpolation k-mesh"							)
		!~~~~~~~~~~~~
		!
		![Fermi]
		call CFG_add(my_cfg,	"Fermi%N_eF"					,	1				,	"number of fermi energys to test"				)
		call CFG_add(my_cfg,	"Fermi%eF_min"					,	1.36_dp			,	"minimum fermi energy( in eV)"					)
		call CFG_add(my_cfg,	"Fermi%eF_max"					,	1.36_dp			,	"maximum fermi energy( in eV)"					)
		call CFG_add(my_cfg,	"Fermi%Tkelvin"					,	0.0_dp			,	"Temperature"									)				
		call CFG_add(my_cfg,	"Fermi%N_eta_smr"				,	201				,	"size of smearing linspace")
		call CFG_add(my_cfg,	"Fermi%eta_smr_min"				,	0.0_dp			,	"min smearing value in eV")
		call CFG_add(my_cfg,	"Fermi%eta_smr_max"				,	2.0_dp			,	"max smearing value in eV")
		call CFG_add(my_cfg,	"Fermi%kuboTol"					,	1e-3_dp			,	"numerical tolearnce for KUBO formulas"			)
		!~~~~~~~~~~~~
		!
		![mep]
		call CFG_add(my_cfg,	"MEP%valence_bands"				,	1				,	"number of valence_bands"						)
		call CFG_add(my_cfg,	"MEP%do_write_mep_bands"		,	.True.			,	"write mep tensor band resolved"				)
		!~~~~~~~~~~~~
		!
		![Laser]
		call CFG_add(my_cfg,	"Laser%N_hw"					,	1				,	"points to probe in interval"					)
		call CFG_add(my_cfg,	"Laser%hw_min"					,	1.55_dp			,	"min energy of incoming light"					)
		call CFG_add(my_cfg,	"Laser%hw_max"					,	1.55_dp			,	"max energy of incoming light"					)
		!~~~~~~~~~~~~		
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		write(*,*)	"[write_input_cfg_file]: BEGIN input file:"
		call CFG_write(my_cfg, "stdout", hide_unused=.true.) ! Write to stdout
		write(*,*)	"[write_input_cfg_file]: END input file."
		call CFG_write(my_cfg,fname_cfg)
		write(*,*)	"[write_input_cfg_file]: wrote input to file "//fname_cfg
		!
		!
		return
	end subroutine




















end program