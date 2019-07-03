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
										get_cart_kpt,							&
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
		module procedure z4_tens_normalize_k_int
		module procedure z5_tens_normalize_k_int
		module procedure z6_tens_normalize_k_int
	end interface normalize_k_int


	integer			::		mp_grid(3)
	real(dp)		::		recip_latt(3,3), bz_vol, unit_vol
	logical			::		recip_latt_set = .false.

contains




	subroutine print_kSpace_info()
		real(dp)		:: dV
		
		if( recip_latt_set) then
			! debug BZ volume 
			dV		=	bz_vol - 8.0_dp * pi_dp**3 / unit_vol
			if(	abs(dV)	> 1e-1_dp) then
				write(*,*)	"[set_recip_latt]: 	WARNING	KSPACE CORRUPTED!! bz_vol /= (8 pi**3/ unit_vol)"
			else
				write(*,'(a,i7.7,a,i5,a,i5,a,i5,a,a,f8.4,a,f8.4,a)')	&
									"[#",mpi_id,";k_space]: k-space setup done (mp_grid=",&
									mp_grid(1),"x",mp_grid(2),"x",mp_grid(3),"). ",&
									"unit_vol=,",unit_vol,",(a0**3) -> bz_vol=",bz_vol," (1/a0**3)"
			end if
		else
			write(*,'(a,i7.7,a)')	&
				"[#",mpi_id,";k_space]: ERROR reciprocal lattice was not initalized. Pls contact dev.."
		end if
		!
		!
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
		call check_lattice()
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
			call check_lattice()
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
			call check_lattice()
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
			call check_lattice()
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
			call check_lattice()
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
		call check_lattice()
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
			call check_lattice()
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
			call check_lattice()
			n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
			if(n_k > 0	)		tens		=	tens		 /	abs(	unit_vol	* real(n_k,dp)	)
		end if
		!
		return
	end subroutine
!------
	subroutine z4_tens_normalize_k_int(tens)
		complex(dp),	allocatable,	intent(inout)	::	tens(:,:,:,:)
		integer							:: 	n_k
		!
		if(allocated(tens))	then
			call check_lattice()
			n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
			if(n_k > 0	)		tens		=	tens		 /	abs(	unit_vol	* real(n_k,dp)	)
		end if
		!
		return
	end subroutine
!------
subroutine z5_tens_normalize_k_int(tens)
	complex(dp),	allocatable,	intent(inout)	::	tens(:,:,:,:,:)
	integer							:: 	n_k
	!
	if(allocated(tens))	then
		call check_lattice()
		n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
		if(n_k > 0	)		tens		=	tens		 /	abs(	unit_vol	* real(n_k,dp)	)
	end if
	!
	return
end subroutine
!------
subroutine z6_tens_normalize_k_int(tens)
	complex(dp),	allocatable,	intent(inout)	::	tens(:,:,:,:,:,:)
	integer							:: 	n_k
	!
	if(allocated(tens))	then
		call check_lattice()
		n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
		if(n_k > 0	)		tens		=	tens		 /	abs(	unit_vol	* real(n_k,dp)	)
		
	end if
	!
	return
end subroutine




subroutine check_lattice()
	if (.not. recip_latt_set) &
		write(*,*)	"[tens_normalize_k_int]: WARNING the reciprocal lattice was never set! Tensor will be normalized with volume 1"
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
			write(*,'(a,i7.7,a)')		"[#",mpi_id,"; kspace_allocator]: allocated position operator (will use Berry connenction & curv)"	
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
		integer			::	i, j
		!	
		!	see eq.(2)		PRB, 13, 5188 (1976)
		!
		real(dp), intent(in)	::	a_latt(3,3)
		real(dp)				::	a1(3), a2(3), a3(3),	&
									b1(3), b2(3), b3(3)
		!
		recip_latt_set	=	.false.
		!
		a1(:)			=	a_latt(:,1)
		a2(:)			=	a_latt(:,2)
		a3(:)			=	a_latt(:,3)
		!get Unit cell volume
		unit_vol		=	dot_product(	crossP( a1(:) , a2(:)	)		,	a3(:)	)	 
		!
		!setup reciprocal lattice vectors
		b1(:)	= 2.0_dp * pi_dp * crossP( a2(:) , a3(:) ) / unit_vol
		b2(:)	= 2.0_dp * pi_dp * crossP( a3(:) , a1(:) ) / unit_vol
		b3(:)	= 2.0_dp * pi_dp * crossP( a1(:) , a2(:) ) / unit_vol
		!
		! BZ volume
		bz_vol	=	dot_product(		crossP( b1(:), b2(:))	, b3(:)		)
		!
		!	CPY TO TARGET 
		recip_latt(:,1)	=	b1(:)
		recip_latt(:,2)	=	b2(:)
		recip_latt(:,3)	=	b3(:)
		recip_latt_set	=	.true.
		! -----
		!
		!	DEBUG
		do i = 1, 3 
			if( 	abs(dot_product(recip_latt(i,:),a_latt(i,:))-2.0_dp * pi_dp)		> 1e-3_dp	) &
				write(*,'(a,i7.7,a,i1,a,i1,a,f8.4,a)')	"[#",mpi_id,&
							"set_recip_latt]:	WARNING b",i,".a",i," =",dot_product(recip_latt(i,:),a_latt(i,:)),'/= 2pi'
			!
			!
			do j = 1, 3
				if(i==j) cycle
				if( 	abs(dot_product(recip_latt(i,:),a_latt(j,:)))		> 1e-3_dp	) &
					write(*,'(a,i7.7,a,i1,a,i1,a,f8.4,a)')	"[",mpi_id,&
							"set_recip_latt]:	WARNING b",i,".a",j," =",dot_product(recip_latt(i,:),a_latt(j,:)),'/= 0'	
			end do
			!
			!
		end do
		!
		if(mpi_id	==	0	)	&
			call print_recip_latt()
		!
		return
	end subroutine

	subroutine print_recip_latt()
		if(mpi_id	==	0	)	&
			write(*,*)	"reciprocal lattice:"
			write(*,*)	recip_latt

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


	function get_rel_kpt(qi_idx, qix,qiy,qiz, kmax) result(kpt)
		integer,				intent(in)	::	qi_idx, qix, qiy, qiz
		real(dp),	optional,	intent(in)	::	kmax
		real(dp),		allocatable			::	kpt(:)
		real(dp)							::	kmax_int = 1.0_dp
		logical,		save				::	k_cut_msg=.False.
		!
		if(check_kpt_idx(qi_idx, qix,qiy,qiz))	then
			allocate(kpt(3))
			!
			if(present(kmax))	&
				kmax_int	=	kmax
			!
			if (.not. k_cut_msg .and. abs(kmax_int-1.0_dp)>1e-3_dp)	&
				write(*,'(a,i5.5,a,e10.3,a)')	"[#",mpi_id,";get_rel_kpt]:	WARNING	k_cutoff=",kmax_int,&
								"	(if you dont know what this means you should fuck off)"
			k_cut_msg	=	.true.
			!
			if( mp_grid(1) > 0 .and. mp_grid(2) > 0 .and. mp_grid(3) > 0 ) then
				!	mp mesh
				kpt(1)	=	kmax_int *	(	 2.0_dp*real(qix,dp)	- real(mp_grid(1),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(1),dp) )
				kpt(2)	=	kmax_int *	(	 2.0_dp*real(qiy,dp)	- real(mp_grid(2),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(2),dp) )
				kpt(3)	=	kmax_int *	(	 2.0_dp*real(qiz,dp)	- real(mp_grid(3),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(3),dp) )
				!
				!	Wx mesh:
				!kpt(1)	=	real(qix-1,dp) / real(mp_grid(1))
				!kpt(2)	=	real(qiy-1,dp) / real(mp_grid(2))
				!kpt(3)	=	real(qiz-1,dp) / real(mp_grid(3))
			end if
		else 
			stop "[get_rel_kpt]: ERROR got illegal kpt idx"
		end if
		!
		!
		return
	end function


	function get_cart_kpt(a_latt, kpt) result(kpt_cart)
		integer			::	j
		!
		!	use the a_latt instead of recip_latt
		!
		real(dp),		intent(in)	::	a_latt(3,3), kpt(3)
		real(dp)					::	kpt_cart(3)
		!
		if(	recip_latt_set) then
			kpt_cart	=	matmul(	recip_latt, 	kpt	)
			!
		else
			write(*,*)	"[kspace/get_cart_kpt]:	ERROR get_cart_kpt was called without lattice beeing set"
			STOP
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