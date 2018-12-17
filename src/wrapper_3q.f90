module wrapper_3q
	use constants,			only:		sp, dp, aUtoEv
	use mpi_comm,			only:		mpi_id
	use m_setham_FeMn,		only:		read_inp,	init_ham
	use input_paras,		only:		debug_mode

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
		!
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^				
		!			PRINT HAM																 |
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		!
		!	sp
		if(debug_mode)	then
			write(*,'(a,i3,a)')	'[#',mpi_id,';	get_ham]	finished hamiltonian setup H_sp:'
			do row = 1, num_wann
				write(*,"(100(a,f5.2,a,f5.2,a))")	(	"(",real(H_sp(row,clm)),"+i* ",imag(H_sp(row,clm)),")	"	,	clm = 1, num_wann		)
			end do
			!
			write(*,'(a,i3,a)')	'[#',mpi_id,';	get_ham]	got velocities:'
			write(*,*)	"	n	| 	m 	||		vx 			|			vy 		|			vz				"
			write(*,*)	"-------------------------------------------------------------------------------------------------"
			do m =	1, size(V_dp,3)
				do n = 1, size(V_dp,2)
				!	
				write(*,'(a,i4,a,i4,a)',advance="no")			" ",n," | ",m, "	||	"
				write(*,'(a,f7.3,a,f7.3,a)',advance="no")		" (",dreal(	V_dp(1,n,m)	),"+i*",aimag( V_dp(1,n,m)),")	| "
				write(*,'(a,f7.3,a,f7.3,a)',advance="no")		" (",dreal(	V_dp(2,n,m)	),"+i*",aimag( V_dp(2,n,m)),")	| "
				write(*,'(a,f7.3,a,f7.3,a)')					" (",dreal(	V_dp(3,n,m)	),"+i*",aimag( V_dp(3,n,m)),")	| "

				!
				end do
				write(*,*)	"------"
			end do
			write(*,*)	"-------------------------------------------------------------------------------------------------"
		end if
		!
		!	dp
		!write(*,'(a,i3,a)')	'[#',mpi_id,';	get_ham]	finished hamiltonian setup H_dp:'
		!do row = 1, num_wann
		!	write(*,"(100(a,f5.2,a,f5.2,a))")	(	"(",dreal(H_dp(row,clm)),"+i* ",aimag(H_dp(row,clm)),")	"	,	clm = 1, num_wann		)
		!end do
		!
		!
	end subroutine











end module wrapper_3q



