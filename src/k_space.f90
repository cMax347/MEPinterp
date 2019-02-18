module k_space
	use constants, 			only:		dp, pi_dp
	use mpi_community,		only:		mpi_id
	use matrix_math,		only:		crossP

	implicit none

	public						::		set_mp_grid,		get_mp_grid,		&
										set_recip_latt, 	get_recip_latt,		&
										get_kpt_idx,							&
										check_kpt_idx,							&
										get_rel_kpt,							&		
										get_bz_vol,								&
										kspace_allocator,						&
										normalize_k_int,						&
										print_kSpace_info
	private	

	save

	interface normalize_k_int
		module procedure d0_tens_normalize_k_int
		module procedure d2_tens_normalize_k_int
		module procedure d3_tens_normalize_k_int
		module procedure d4_tens_normalize_k_int
		module procedure d5_tens_normalize_k_int
		!
		module procedure z0_tens_normalize_k_int
		module procedure z2_tens_normalize_k_int
		module procedure z3_tens_normalize_k_int
	end interface normalize_k_int


	integer			::		mp_grid(3)
	real(dp)		::		recip_latt(3,3), bz_vol, unit_vol

contains




	subroutine print_kSpace_info()
		write(*,'(a,i3,a,i4,a,i4,a,i4,a,a,f8.4,a)')	&
							"[#",mpi_id,";k_space]: k-space setup done (mp_grid=",&
							mp_grid(1),"x",mp_grid(2),"x",mp_grid(3),"). ",&
							"bz_vol=",bz_vol," (1/a0)"

		return
	end subroutine

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------		
!
!			INTERFACE:		normalize_k_int
!				
!			Vbz	=	(2pi)^3/ Vcell			
!
!
!			int 	Vcell / (2pi)^3		d^3k  		---->		1 /( N_k)	sum_k
!
!-----------------------------------------------------------------------------------
	subroutine d0_tens_normalize_k_int(scalar)
		real(dp),	intent(inout)	::	scalar
		integer						::	n_k
		!
		n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
		!
		!
		if(n_k > 0	)		scalar	= 	scalar	 /	abs(	unit_vol	* real(n_k,dp)	)
		
		!
		return
	end subroutine
!------
	subroutine d2_tens_normalize_k_int(tens)
		real(dp),	allocatable,	intent(inout)	::	tens(:,:)
		integer						:: 	n_k
		!
		if(allocated(tens)) then
			n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
			if(n_k > 0	)		tens		=	tens		 /	abs(	unit_vol	* real(n_k,dp)	)
		end if
		!
		return
	end subroutine
!------
	subroutine d3_tens_normalize_k_int(tens)	
		real(dp),	allocatable,	intent(inout)	::	tens(:,:,:)
		integer						:: 	n_k
		!
		if(allocated(tens)) then
			n_k	=	mp_grid(1) * mp_grid(2) * mp_grid(3)
			if(n_k > 0	)		tens	=	tens	 /	abs(	unit_vol	* real(n_k,dp)	)
		end if
		!
		return
	end subroutine
!------
	subroutine d4_tens_normalize_k_int(tens)	
		real(dp),	allocatable,	intent(inout)	::	tens(:,:,:,:)
		integer						:: 	n_k
		!
		if(allocated(tens)) then
			n_k	=	mp_grid(1) * mp_grid(2) * mp_grid(3)
			if(n_k > 0	)		tens	=	tens	 /	abs(	unit_vol	* real(n_k,dp)	)
		end if
		!
		return
	end subroutine
!------
	subroutine d5_tens_normalize_k_int(tens)	
		real(dp),	allocatable,	intent(inout)	::	tens(:,:,:,:,:)
		integer						:: 	n_k
		!
		if(allocated(tens)) then
			n_k	=	mp_grid(1) * mp_grid(2) * mp_grid(3)
			if(n_k > 0	)		tens	=	tens	 /	abs(	unit_vol	* real(n_k,dp)	)
		end if
		!
		return
	end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~
!------
!^^^^^^^^^^^^^^^^^^^^^^^^^
	subroutine z0_tens_normalize_k_int(scalar)
		complex(dp),	intent(inout)	::	scalar
		integer							::	n_k
		!
		n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
		!
		!
		if(n_k > 0	)		scalar	= 	scalar	 /	abs(	unit_vol	* real(n_k,dp)	)
		!
		return
	end subroutine
!------
	subroutine z2_tens_normalize_k_int(tens)
		complex(dp),	allocatable,	intent(inout)	::	tens(:,:)
		integer							:: 	n_k
		!
		if(allocated(tens))	then
			n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
			if(n_k > 0	)		tens		=	tens		 /	abs(	unit_vol	* real(n_k,dp)	)
		end if
		!
		return
	end subroutine
!------
	subroutine z3_tens_normalize_k_int(tens)
		complex(dp),	allocatable,	intent(inout)	::	tens(:,:,:)
		integer							:: 	n_k
		!
		if(allocated(tens))	then
			n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
			if(n_k > 0	)		tens		=	tens		 /	abs(	unit_vol	* real(n_k,dp)	)
		end if
		!
		return
	end subroutine
	












!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

	subroutine kspace_allocator(H_tb, r_tb, en_k, V_ka, A_ka, Om_kab)
		!	allocate the k-space arrays, depending on what is allocated in real space
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
			write(*,'(a,i3,a)')		"[#",mpi_id,"; kspace_allocator]: allocated position operator (will use Berry connenction & curv)"	
		end if
		!
		return
	end subroutine





!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

	subroutine set_mp_grid(input_grid)
		!	set the monkhorst pack grid
		integer,		intent(in)		::	input_grid(3)
		mp_grid	=	input_grid
		return
	end subroutine

	function get_mp_grid() result(this_mp_grid)
		integer			::	this_mp_grid(3)
		this_mp_grid	=	mp_grid
		return
	end function

	




	subroutine set_recip_latt(a_latt)
		!	
		!	see eq.(2)		PRB, 13, 5188 (1976)
		!
		real(dp), intent(in)	::	a_latt(3,3)
		real(dp)				::	a1(3), a2(3), a3(3),	&
									b1(3), b2(3), b3(3)
		!
		!
		a1(1:3)		=	a_latt(1,1:3)
		a2(1:3)		=	a_latt(2,1:3)
		a3(1:3)		=	a_latt(3,1:3)
		!get Unit cell volume
		unit_vol	=	dot_product(	crossP( a1(1:3) , a2(1:3)	)		,	a3(1:3)	)	 
		!
		!setup reciprocal lattice vectors
		b1(1:3)	= 2.0_dp * pi_dp * crossP( a2(1:3) , a3(1:3) ) / unit_vol
		b2(1:3)	= 2.0_dp * pi_dp * crossP( a3(1:3) , a1(1:3) ) / unit_vol
		b3(1:3)	= 2.0_dp * pi_dp * crossP( a1(1:3) , a2(1:3) ) / unit_vol
		!
		! BZ volume
		bz_vol	=	dot_product(		crossP( b1(:), b2(:))	, b3(:)		)
		!
		!	CPY TO TARGET 
		recip_latt(1,:)	=	b1(:)
		recip_latt(2,:)	=	b2(:)
		recip_latt(3,:)	=	b3(:)
		!
		return
	end subroutine



	function get_recip_latt() result( this_recip_latt )
		real(dp)		::	this_recip_latt(3,3)
		this_recip_latt	= recip_latt
		return
	end function

	function get_bz_vol() result(this_bz_vol)
		real(dp)		::	this_bz_vol
		this_bz_vol=bz_vol
		return 
	end function



	integer pure function get_kpt_idx(qix, qiy, qiz)
		integer,	intent(in)	::	qix, qiy, qiz
		!
		!	use start at 0 convention and ad 1 to switch to array starts at 1 convention
		get_kpt_idx	=	1	+	(qix-1) + mp_grid(2)	* ( (qiy-1) + mp_grid(3) * (qiz-1)	)
		!
		return
	end function


	function get_rel_kpt(qi_idx, qix,qiy,qiz) result(kpt)
		integer,	intent(in)	::	qi_idx, qix, qiy, qiz
		real(dp),	allocatable	::	kpt(:)
		!
		if(check_kpt_idx(qi_idx, qix,qiy,qiz))	then
			allocate(kpt(3))
			!
			if( mp_grid(1) > 0 .and. mp_grid(2) > 0 .and. mp_grid(3) > 0 ) then
				kpt(1)	=	(	 2.0_dp*real(qix,dp)	- real(mp_grid(1),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(1),dp) )
				kpt(2)	=	(	 2.0_dp*real(qiy,dp)	- real(mp_grid(2),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(2),dp) )
				kpt(3)	=	(	 2.0_dp*real(qiz,dp)	- real(mp_grid(3),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(3),dp) )
			end if
		else 
			stop "[get_rel_kpt]: ERROR got illegal kpt idx"
		end if
		!
		!
		return
	end function


	logical function check_kpt_idx(	qi_idx	,qix,qiy,qiz)
		integer,	intent(in)	::	qi_idx, qix, qiy, qiz
		!
		!	debug check	
		if	(					qi_idx 		< 				0 						&
					.or.		qi_idx		> 	mp_grid(1)*mp_grid(2)*mp_grid(3)	) 	then
			!
			stop "[check_kpt_idx]: ERROR - idx out of bounds"
		end if
		!
		check_kpt_idx	=	(	qi_idx	-	get_kpt_idx(qix,qiy,qiz) )	==	0
		!
		return
	end function



end module k_space