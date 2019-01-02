module wrapper_3q
	use constants,			only:		sp, dp, aUtoEv
	use mpi_comm,			only:		mpi_id
	use m_setham_FeMn,		only:		read_inp,	init_ham

	implicit none


	private
	public					::		get_ham

	interface rm_num_noise
        module procedure   d_zet_zero
        module procedure   z_zet_zero
    end interface rm_num_noise


	contains



	!PUBLIC:
	subroutine get_ham(rel_kpt_dp,	H_dp, V_dp)
		real(dp),		intent(in)							::	rel_kpt_dp(3)
		complex(dp),	allocatable,	intent(inout)		::	H_dp(:,:),	V_dp(:,:,:)
		complex(sp),	allocatable							::	H_sp(:,:),	vx_sp(:,:), vy_sp(:,:),	vz_sp(:,:), &
																conn_dummy(:,:,:), curv_dummy(:,:,:)
		integer												::	num_wann, row, clm, n, m, x
		real(dp)											::	single_prec
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
		H_dp				=	 H_sp(:,:) 	/	aUtoEv
		V_dp(	1	,:,:)	=	vx_sp(:,:)	/	aUtoEv
		V_dp(	2	,:,:)	=	vy_sp(:,:)	/	aUtoEv
		V_dp(	3	,:,:)	=	vz_sp(:,:)	/	aUtoEv
		!
		!	
		!	make Ham hermitian
		do m = 1, size(H_dp,2)
			do n = 1, size(H_dp,1)
				if( n > m ) 	H_dp(n,m)	=	conjg(	H_dp(m,n)	)	
			end do
		end do
		!
		!	try to avoid 1e-9 hoppings in the hamiltonian (allows for better comparability)
		single_prec		=	1e-8_dp
		do m = 1, size(H_dp,2)
			do n = 1, size(H_dp,1)
				call rm_num_noise(		H_dp(	n,m)	)
				!
				do x = 1, 3
					call rm_num_noise(	V_dp(x,	n,m)	)
				end do
			end do
		end do
		!
		!
	end subroutine









	pure subroutine z_zet_zero(	z)
		complex(dp),	intent(inout)		::	z
		real(dp)							::	re_z,	im_z	
		!
		re_z	=	real(	z	,dp	)
		im_z	=	aimag(	z		)
		!
		call d_zet_zero(	re_z	)
		call d_zet_zero(	im_z	)
		!
		z	=	cmplx(re_z, im_z, dp)
		!
	end subroutine


	pure subroutine d_zet_zero(	d)
		real(dp),	intent(inout)			::	d
		!
		if(	abs(d)	< 1e-8_dp)		d	=	0.0_dp
	end subroutine








end module wrapper_3q



