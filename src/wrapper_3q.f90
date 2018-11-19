module wrapper_3q
	use constants,			only:		sp, dp
	use mpi_comm,			only:		mpi_id
	use m_setham_FeMn,		only:		read_inp,	init_ham



	implicit none


	private
	public					::		get_ham


	contains



	!PUBLIC:
	subroutine get_ham(kpt_dp,	H_dp, V_dp)
		real(dp),		intent(in)							::	kpt_dp(3)
		complex(dp),	allocatable,	intent(inout)		::	H_dp(:,:),	V_dp(:,:,:)
		complex(sp),	allocatable							::	H_sp(:,:),	vx_sp(:,:), vy_sp(:,:),	vz_sp(:,:)
		integer												::	num_wann	
		!
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			ALLOCATE																  |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		num_wann	=	8
		allocate(		H_sp(	num_wann,	num_wann)		)
		allocate(		H_dp(	num_wann,	num_wann)		)
		!
		if(	allocated(	V_dp)	) then
			allocate(	vx_sp(	num_wann,	num_wann ))				
			allocate(	vy_sp(	num_wann,	num_wann ))
			allocate(	vz_sp(	num_wann,	num_wann ))
		end if
		!
		!
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
		write(*,*)	'[#',mpi_id,';	get_ham]	Init 3q state k-space ham for kpt=',kpt_dp		
		if(	allocated(V_dp)	) 	then
			!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
			!			ALSO SETUP VELOCITIES												  |
			!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			call init_ham(	real(kpt_dp,sp), 	num_wann,	H_sp, 	vx_sp, vy_sp, vz_sp	)
		else
			!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
			!			SETUP HAM ONLY														  |
			!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			call init_ham( 	real(kpt_dp,sp), 	num_wann,	H_sp						)
		end if
		!
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			CONVERT SINGLE TO DOUBLE												 |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		H_dp	=	H_sp
		if(	allocated(V_dp))	then
			V_dp(	1	,:,:)	=	vx_sp(:,:)
			V_dp(	2	,:,:)	=	vy_sp(:,:)
			V_dp(	3	,:,:)	=	vz_sp(:,:)	
		end if
		!
		!
		!
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			PRINT HAM																 |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		write(*,'(a,i3,a)')	'[#',mpi_id,';	get_ham]	finished hamiltonian setup H_sp:'
		write(*,*)	H_sp
		write(*,'(a,i3,a)')	'[#',mpi_id,';	get_ham]	finished hamiltonian setup H_dp:'
		write(*,*)	H_dp
	end subroutine











end module wrapper_3q



