module wrapper_3q
	use constants,			only:		sp, dp, aUtoEv
	use mpi_comm,			only:		mpi_id
	use m_setham_FeMn,		only:		read_inp,	init_ham

	implicit none


	private
	public					::		get_ham


	contains



	!PUBLIC:
	subroutine get_ham(rel_kpt_dp,	H_dp, V_dp)
		real(dp),		intent(in)							::	rel_kpt_dp(3)
		complex(dp),	allocatable,	intent(inout)		::	H_dp(:,:),	V_dp(:,:,:)
		complex(sp),	allocatable							::	H_sp(:,:),	vx_sp(:,:), vy_sp(:,:),	vz_sp(:,:), &
																conn_dummy(:,:,:), curv_dummy(:,:,:)
		integer												::	num_wann, row, clm, n, m	
		!
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			ALLOCATE																  |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		num_wann	=	8
		allocate(		H_sp(		num_wann,	num_wann)		)
		if(.not. allocated(H_dp))	allocate(		H_dp(		num_wann,	num_wann)		)
		if(.not. allocated(V_dp))	allocate(		V_dp(	3,	num_wann,	num_wann)		)
		!
		allocate(	vx_sp(	num_wann,	num_wann ))				
		allocate(	vy_sp(	num_wann,	num_wann ))
		allocate(	vz_sp(	num_wann,	num_wann ))
		!
		!
		H_sp	=	0.0
		vx_sp	=	0.0
		vy_sp	=	0.0
		vz_sp	=	0.0



		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			CONFIGURE HOPPINGS, ETC.				   								 |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		write(*,'(a,i3,a)')	'[#',mpi_id,';	get_ham]	try to read 3q_state input file...'
		call read_inp(mpi_id, .true.)
		!
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			HAM SETUP																  |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		write(*,'(a,i3,a,f7.3,a,f7.3,a,f7.3,a)')	'[#',mpi_id,													&
													';	get_ham]	Init 3q state k-space ham for rel. kpt= (',		&
													rel_kpt_dp(1),", ",rel_kpt_dp(2),", ",rel_kpt_dp(3)	, ")."

		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			ALSO SETUP VELOCITIES												  |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		call init_ham(	real(rel_kpt_dp,sp), 	num_wann,	H_sp, 	vx_sp, vy_sp, vz_sp	)
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			CONVERT SINGLE TO DOUBLE												 |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		!		jan seems to use atomic units, however energies are in (eV) ?!
		!
		H_dp				=	real(	 H_sp(:,:)	,dp) 	/	aUtoEv
		V_dp(	1	,:,:)	=	real(	vx_sp(:,:)	,dp)	/	aUtoEv
		V_dp(	2	,:,:)	=	real(	vy_sp(:,:)	,dp)	/	aUtoEv
		V_dp(	3	,:,:)	=	real(	vz_sp(:,:)	,dp)	/	aUtoEv
		!
		!
	end subroutine











end module wrapper_3q



