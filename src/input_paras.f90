module input_paras
	use m_config
#ifdef USE_MPI
	use mpi
#endif
	use matrix_math,				only:		crossP
	use constants,					only:		dp, fp_acc, pi_dp, aUtoEv
	use mpi_comm,					only:		mpi_id, mpi_root_id, mpi_nProcs, ierr
	use k_space,					only:		set_recip_latt, set_mp_grid




	implicit none

	private
	public								::		&
												!routines		
												init_parameters, 	my_mkdir,						 		 		&
												!dirs
												w90_dir, 															&
												raw_dir,															&
												out_dir,															&
												eig_out_dir,														&
												velo_out_dir,														&
												mep_out_dir, ahc_out_dir, 	opt_out_dir,	gyro_out_dir,			&
												!jobs
												plot_bands,															&
												use_mpi,															&
												debug_mode,															&
												do_gauge_trafo,														&
												do_write_velo,														&
												use_R_float,														&
												do_write_mep_bands,													&
												do_mep, do_ahc, do_kubo, do_opt, do_gyro,							&
												!vars
												seed_name,	valence_bands,											&
												a_latt, kubo_tol, unit_vol,											&
												hw, phi_laser, eFermi, T_kelvin, i_eta_smr



	
	save

	
	integer						::	valence_bands
	character(len=3)			:: 	seed_name
	character(len=9)			::	w90_dir	="w90files/"
	character(len=4)			::	raw_dir ="raw/"
	character(len=4)			::	out_dir	="out/"
	character(len=10)			:: 	velo_out_dir
	character(len=11)			::	eig_out_dir
	character(len=9)			::	mep_out_dir	
	character(len=9)			::	ahc_out_dir
	character(len=9)			::	opt_out_dir
	character(len=10)			::	gyro_out_dir	
	logical						::	plot_bands, do_gauge_trafo, 	&
									use_R_float,					&
									do_write_velo,					&
									do_write_mep_bands,				&
									debug_mode,	use_mpi,			&
									do_mep, do_ahc, do_kubo, do_opt, do_gyro
	real(dp)					::	a_latt(3,3), a0, unit_vol,		&
									kubo_tol,						&
									hw, eFermi, T_kelvin			
	complex(dp)					::	i_eta_smr, phi_laser




	contains




!public
	logical function init_parameters()
		type(CFG_t) 			:: 	my_cfg
		real(dp)				::	a1(3), a2(3), a3(3), eta, laser_phase
		integer					::	mp_grid(3)
		logical					::	input_exist
		!
		use_mpi	= .false.
#ifdef USE_MPI
		use_mpi = .true.
#endif
		if( 	.not. use_mpi 		.and.		 mpi_id	/= 0			)	then
			 write(*,'(a,i3,a)')		'[#',mpi_id,';init_parameters]:	hello, I am an unexpected MPI thread !!!1!1!!!1!!1 '
		end if
		!
		velo_out_dir	=	out_dir//"/velo/"
		mep_out_dir 	=	out_dir//"/mep/"	
		ahc_out_dir		=	out_dir//"/ahc/"
		opt_out_dir 	=	out_dir//"/opt/"	
		gyro_out_dir	=	out_dir//"/gyro/"
		eig_out_dir		=	raw_dir//"/eigen/"
		!ROOT READ
		if(mpi_id == mpi_root_id) then
			inquire(file="./input.cfg",exist=input_exist)
			!
			if( input_exist)	then		
				!OPEN FILE
				call CFG_read_file(my_cfg,"./input.cfg")
				!
				![methods]
				call CFG_add_get(my_cfg,	"jobs%plot_bands"				,	plot_bands			,	"if true do a bandstructure run"	)
				call CFG_add_get(my_cfg,	"jobs%debug_mode"				,	debug_mode			,	"switch aditional debug tests in code")
				call CFG_add_get(my_cfg,	"jobs%R_vect_float"				,	use_R_float			,	"the R_cell vector is now real (else: integer)")
				call CFG_add_get(my_cfg,	"jobs%do_write_velo"			,	do_write_velo		,	"write formatted velocity files at each kpt")
				call CFG_add_get(my_cfg,	"jobs%do_mep"					,	do_mep				,	"switch response tensor calc")
				call CFG_add_get(my_cfg,	"jobs%do_kubo"					,	do_kubo				,	"switch response tensor calc")
				call CFG_add_get(my_cfg,	"jobs%do_ahc"					,	do_ahc				,	"switch response tensor calc")
				call CFG_add_get(my_cfg,	"jobs%do_opt"					,	do_opt				,	"switch response tensor calc")
				call CFG_add_get(my_cfg,	"jobs%do_gyro"					,	do_gyro				,	"switch response tensor calc")



				![unitCell]
				call CFG_add_get(my_cfg,	"unitCell%a1"      				,	a1(1:3)  	  		,	"a_x lattice vector"				)
				call CFG_add_get(my_cfg,	"unitCell%a2"      				,	a2(1:3)  	  		,	"a_y lattice vector"				)
				call CFG_add_get(my_cfg,	"unitCell%a3"      				,	a3(1:3)  	  		,	"a_z lattice vector"				)
				call CFG_add_get(my_cfg,	"unitCell%a0"					,	a0					,	"lattice scaling factor "			)
				!
				a_latt(1,1:3)	= a1(1:3)
				a_latt(2,1:3)	= a2(1:3)
				a_latt(3,1:3)	= a3(1:3)
				a_latt			=	a0 * a_latt
				!
				!
				![wannInterp]
				call CFG_add_get(my_cfg,	"wannInterp%doGaugeTrafo"		,	do_gauge_trafo		,	"switch (W)->(H) gauge trafo"		)
				call CFG_add_get(my_cfg,	"wannInterp%mp_grid"			,	mp_grid(1:3)		,	"interpolation k-mesh"				)
				call CFG_add_get(my_cfg,	"wannInterp%seed_name"			,	seed_name			,	"seed name of the TB files"			)
				!
				![mep]
				call CFG_add_get(my_cfg,	"MEP%valence_bands"				,	valence_bands		,	"number of valence_bands"			)
				call CFG_add_get(my_cfg,	"MEP%do_write_mep_bands"		,	do_write_mep_bands	,	"write mep tensor band resolved"	)
				![Fermi]
				call CFG_add_get(my_cfg,	"Kubo%kuboTol"					,	kubo_tol			,	"numerical tolearnce for KUBO"		)
				call CFG_add_get(my_cfg,	"Kubo%hw"						,	hw					,	"energy of incoming light"			)
				call CFG_add_get(my_cfg,	"Kubo%laser_phase"				,	laser_phase			,	"euler angle of phase shift of driving E-field")
				call CFG_add_get(my_cfg,	"Kubo%eFermi"					,	eFermi				,	"set the Fermi energy"				)
				call CFG_add_get(my_cfg,	"Kubo%Tkelvin"					,	T_kelvin			,	"Temperature"						)				
				call CFG_add_get(my_cfg,	"Kubo%eta_smearing"				,	eta					,	"smearing for optical conductivty"	)
				!
				! 	unit conversion
				hw			=	hw 		/ 	aUtoEv
				eFermi		= 	eFermi	/	aUtoEv
				eta			=	eta		/	aUtoEv
				!
				!	derived constants
				i_eta_smr	=	cmplx(0.0_dp,	eta	,dp)
				phi_laser	=	cmplx(cos(pi_dp*laser_phase),sin(pi_dp*laser_phase),dp)
				!
				!
				write(*,*)					""
				write(*,*)					"**********************init_parameter interpretation"	&
											//" (3q state mode most paras will be ignored)******************************"
				write(*,*)					"parallelization with ",mpi_nProcs," MPI threads"
				write(*,'(a,i3,a)')			"[#",mpi_id,";init_parameters]: input interpretation:"
				write(*,*)					"[methods]"
				write(*,*)					"	plot_bands=",plot_bands
				write(*,*)					"	debug_mode=",debug_mode	
				write(*,*)					"	do_write_velo=",do_write_velo	
				write(*,*)					"[unitCell] # a_0 (Bohr radii)	"
				write(*,*)					"	a1=",a1(1:3)
				write(*,*)					"	a2=",a2(1:3)
				write(*,*)					"	a3=",a3(1:3)
				write(*,*)					"	a0=",a0
				write(*,*)					"[wannInterp]"
				write(*,*)					"	seed_name=",seed_name
				write(*,*)					"[mep]"
				write(*,'(a,i4)')			"	val bands=",valence_bands
				write(*,*)					"	do_write_mep_bands=",do_write_mep_bands
				write(*,*)					"[Kubo] # E_h (Hartree)"
				write(*,*)					"	do_gauge_trafo=",do_gauge_trafo
				write(*,*)					"	kuboTol=",kubo_tol
				write(*,*)					"	hw=",hw*aUtoEv," (eV)" 
				write(*,*)					"	laser_phase= pi * ",laser_phase
				write(*,*)					"	phi_laser( exp(i*laser_phase))=",phi_laser
				write(*,*)					"	eFermi=",eFermi*aUtoEv," (eV)"
				write(*,*)					"	T_kelvin=",T_kelvin," (K)"
				write(*,*)					"	eta=",eta*aUtoEv," (eV)"
				write(*,*)					"	i_eta_smr=",i_eta_smr
				write(*,*)					"*********************************************************************************"		
				!
				!make the output folder
				write(*,*)	"*"
				write(*,*)	"----------------------MAKE TARGET DIRECTORIES-------------------------"
				write(*,'(a,i3,a)')			"[#",mpi_id,";init_parameters]: start target mkdir..."
				call my_mkdir(out_dir)
				call my_mkdir(raw_dir)
				if(		 	do_write_velo		)	call my_mkdir(velo_out_dir)	
				if(			debug_mode			)	call my_mkdir(eig_out_dir)
				if(.not. plot_bands) then
					if( do_mep .or. do_kubo		)	call my_mkdir(mep_out_dir)
					if(			do_ahc			)	call my_mkdir(ahc_out_dir)
					if(			do_opt			)	call my_mkdir(opt_out_dir)
					if(			do_gyro			)	call my_mkdir(gyro_out_dir)		
				end if	
				write(*,'(a,i3,a)')			"[#",mpi_id,";init_parameters]: ... directories created"
			else
				write(*,'(a,i3,a)')			"[#",mpi_id,";init_parameters]: could not find input file"
				stop "please provide a input.cfg file"
			end if
			write(*,*)	"*"
			write(*,*)	"----------------------K-SPACE SETUP-------------------------"
			write(*,'(a,i3,a)')			"[#",mpi_id,";init_parameters]: now bcast the input parameters and setup k-space ..."
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
				call MPI_BCAST(		plot_bands		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		debug_mode		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		do_gauge_trafo	,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		do_write_velo	,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_write_mep_bands,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		use_R_float		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_mep 			,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_kubo 		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_ahc 			,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_opt 			,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		do_gyro 		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				!
				call MPI_BCAST(		a_latt			,			9			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		valence_bands	,			1			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		seed_name(:)	,	len(seed_name)		,		MPI_CHARACTER		,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		mp_grid			,			3			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				![KUBO]
				call MPI_BCAST(		kubo_tol		,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD,	ierr)	
				call MPI_BCAST(		hw				,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		eFermi			,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
				call MPI_BCAST(		T_kelvin		,			1			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		i_eta_smr		,			1			,	MPI_DOUBLE_COMPLEX		,		mpi_root_id,	MPI_COMM_WORLD, ierr)
				call MPI_BCAST(		phi_laser		,			1			,	MPI_DOUBLE_COMPLEX		,		mpi_root_id,	MPI_COMM_WORLD, ierr)
			end if
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
		end if
		!
		init_parameters	=	input_exist
		!
		call MPI_BARRIER(MPI_COMM_WORLD, ierr)
		if(mpi_id	== mpi_root_id)		then
			write(*,*)	""
			write(*,'(a,i3,a)')	"[#",mpi_id,";init_parameters]: ...input paramaters broadcasted, k-space setup complete!"
		end if
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
			write(*,'(a,i3,a,a)')	"[#",mpi_id,"; init_parameters]: (fake) created directory ",dir
		!end if
		!
		return
	end subroutine




end module input_paras