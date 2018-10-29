module k_space
	use constants, 			only:		dp, pi_dp, mpi_id
	use matrix_math,		only:		crossP

	implicit none

	public						::		set_mp_grid,		get_mp_grid,		&
										set_recip_latt, 	get_recip_latt,		&
										get_rel_kpt,		get_bz_vol,			&
										normalize_k_int
	private	

	save

	interface normalize_k_int
		module procedure d_scalar_normalize_k_int
		module procedure z_scalar_normalize_k_int
		module procedure d_matrix_normalize_k_int
		module procedure z_matrix_normalize_k_int
	end interface normalize_k_int


	integer			::		mp_grid(3)
	real(dp)		::		recip_latt(3,3), bz_vol

contains

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------		
!
!			INTERFACE:		normalize_k_int
!
!-----------------------------------------------------------------------------------
	subroutine z_matrix_normalize_k_int(mat)
		complex(dp),	intent(inout)	::	mat(:,:)
		integer							:: 	n_k
		!
		n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
		!
		!
		if(n_k > 0	)		mat		=	mat		* bz_vol	/ real(n_k,dp)
		!
		return
	end subroutine
	


	subroutine d_matrix_normalize_k_int(mat)
		real(dp),	intent(inout)	::	mat(:,:)
		integer						:: 	n_k
		!
		n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
		!
		!
		if(n_k > 0	)		mat		=	mat		* bz_vol	/ real(n_k,dp)
		!
		return
	end subroutine

	subroutine d_scalar_normalize_k_int(scalar)
		real(dp),	intent(inout)	::	scalar
		integer						::	n_k
		!
		n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
		!
		!
		if(n_k > 0	)		scalar	= 	scalar	* bz_vol	/ real(n_k,dp)
		
		!
		return
	end subroutine

	subroutine z_scalar_normalize_k_int(scalar)
		complex(dp),	intent(inout)	::	scalar
		integer							::	n_k
		!
		n_k		=	mp_grid(1) * mp_grid(2) * mp_grid(3)
		!
		!
		if(n_k > 0	)		scalar	= 	scalar	* bz_vol	/ real(n_k,dp)
		!
		return
	end subroutine

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

	subroutine set_mp_grid(input_grid)
		integer,		intent(in)		::	input_grid(3)
		mp_grid	=	input_grid
		write(*,'(a,i3,a,(20i5))') 	'[#',mpi_id,';set_mp_grid]: mp_grid set to ', mp_grid
		!						TODO USE MPI_ID HERE AS WELL, ORE REMOVE MESSAGE, ELSE DOES NOT MAKE MUCH SENSE
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
		real(dp)				::	unit_vol,				& 
									a1(3), a2(3), a3(3),	&
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
		write(*,'(a,i3,a,f6.3,a)') 	'[#',mpi_id,';set_recip_latt]: the 1st Brillouin zone volume is	 bz_vol=',bz_vol,'	(1/a_0)^3'
		!
		!	CPY TO TARGET 
		recip_latt(1,:)	=	b1(:)
		recip_latt(2,:)	=	b2(:)
		recip_latt(3,:)	=	b3(:)
		!write(*,*) 				'[set_recip_latt]: recip_latt set to (1/a_0) '
		!write(*,*)				'	',	recip_latt(1,:)
		!write(*,*)				'	',	recip_latt(2,:)
		!write(*,*)				'	',	recip_latt(3,:)
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



	integer function get_rel_kpt(qix, qiy, qiz, kpt)
		integer,	intent(in)	::	qix, qiy, qiz
		real(dp),	intent(out)	::	kpt(3)
		!
		if( mp_grid(1) > 0 .and. mp_grid(2) > 0 .and. mp_grid(3) > 0 ) then
			kpt(1)	=	(	 2.0_dp*real(qix,dp)	- real(mp_grid(1),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(1),dp) )
			kpt(2)	=	(	 2.0_dp*real(qiy,dp)	- real(mp_grid(2),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(2),dp) )
			kpt(3)	=	(	 2.0_dp*real(qiz,dp)	- real(mp_grid(3),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(3),dp) )
		end if
		!
		!	start at 0 convention
		get_rel_kpt	=	(qix-1) + mp_grid(2)	* ( (qiy-1) + mp_grid(3) * (qiz-1)	)
		!
		!	start at 1 convention
		get_rel_kpt	=	get_rel_kpt + 1
		!
		return
	end function



end module k_space