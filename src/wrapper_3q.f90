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
		complex(sp),	allocatable							::	H_sp(:,:),	vx_sp(:,:), vy_sp(:,:),	vz_sp(:,:)
		integer												::	num_wann, row, clm	
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
		write(*,*)	'[#',mpi_id,';	get_ham]	Init 3q state k-space ham for rel. kpt=',rel_kpt_dp		
		if(	allocated(V_dp)	) 	then
			!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
			!			ALSO SETUP VELOCITIES												  |
			!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			call init_ham(	real(rel_kpt_dp,sp), 	num_wann,	H_sp, 	vx_sp, vy_sp, vz_sp	)
		else
			!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
			!			SETUP HAM ONLY														  |
			!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			call init_ham( 	real(rel_kpt_dp,sp), 	num_wann,	H_sp						)
		end if
		!
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			CONVERT SINGLE TO DOUBLE												 |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		H_dp	=	H_sp /	aUtoEv
		V_dp(	1	,:,:)	=	vx_sp(:,:)
		V_dp(	2	,:,:)	=	vy_sp(:,:)
		V_dp(	3	,:,:)	=	vz_sp(:,:)	
		!
		!
		!
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			PRINT HAM																 |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		!
		!	dp
		write(*,'(a,i3,a)')	'[#',mpi_id,';	get_ham]	finished hamiltonian setup H_sp:'
		do row = 1, num_wann
			write(*,"(100(a,f5.2,a,f5.2,a))")	(	"(",real(H_sp(row,clm)),"+i* ",imag(H_sp(row,clm)),")	"	,	clm = 1, num_wann		)
		end do
		!
		!	sp
		write(*,'(a,i3,a)')	'[#',mpi_id,';	get_ham]	finished hamiltonian setup H_dp:'
		do row = 1, num_wann
			write(*,"(100(a,f5.2,a,f5.2,a))")	(	"(",dreal(H_dp(row,clm)),"+i* ",aimag(H_dp(row,clm)),")	"	,	clm = 1, num_wann		)
		end do
		!
		!
	end subroutine











end module wrapper_3q



