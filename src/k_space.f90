module k_space
	use constants, 			only:		dp, pi_dp
	use matrix_math,		only:		crossP

	implicit none

	public						::		set_mp_grid,		get_mp_grid,		&
										set_recip_latt, 	get_recip_latt,		&
										get_rel_kpts,							&
										get_rel_kpt
	private	


	integer			::		mp_grid(3)
	real(dp)		::		recip_latt(3,3)

contains


	subroutine set_mp_grid(input_grid)
		integer,		intent(in)		::	input_grid(3)
		mp_grid	=	input_grid
		write(*,*) '[set_mp_grid]: mp_grid set to ', mp_grid
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
		real(dp)				::	unit_vol, a1(3), a2(3), a3(3)
		!
		!
		a1(1:3)		=	a_latt(1,1:3)
		a2(1:3)		=	a_latt(2,1:3)
		a3(1:3)		=	a_latt(3,1:3)
		!get Unit cell volume
		unit_vol	=	dot_product(	crossP( a1(1:3) , a2(1:3)	)		,	a3(1:3)	)	 
		!
		!
		recip_latt(1,1:3)	= 2.0_dp * pi_dp * crossP( a2(1:3) , a3(1:3) ) / unit_vol
		recip_latt(2,1:3)	= 2.0_dp * pi_dp * crossP( a3(1:3) , a1(1:3) ) / unit_vol
		recip_latt(3,1:3)	= 2.0_dp * pi_dp * crossP( a1(1:3) , a2(1:3) ) / unit_vol
		!
		write(*,*) '[set_recip_latt]: recip_latt set to ', recip_latt
		!
		return
	end subroutine

	function get_recip_latt() result( this_recip_latt )
		real(dp)		::	this_recip_latt(3,3)
		!
		this_recip_latt	= recip_latt
		!
		return
	end function





	subroutine get_rel_kpts(kpt_latt)
		!get relative k-pt following the Monkhorst Pack scheme 
		!	see eq.(4)		PRB, 13, 5188 (1976)
		!
		real(dp),	allocatable, intent(inout)	::	kpt_latt(:,:)
		real(dp)								::	kpt(3)
		integer									::	qix, qiy, qiz, qi_idx, qi_test
		!
		allocate(kpt_latt(3,mp_grid(1)*mp_grid(2)*mp_grid(3)))
		qi_test = 0
		do qiz = 1, mp_grid(3)
			do qiy = 1, mp_grid(2)
				do qix = 1, mp_grid(1)
					qi_idx	= get_rel_kpt(qix,qiy,qiz, kpt	)	
					kpt_latt(:,qi_idx)	=	kpt(:)
					qi_test	= qi_test +1
					if(qi_idx /= qi_test) then
						write(*,'(a,i10,a,i10)')	'[get_rel_kpts]: WARNING k-mesh order! qi_idx=',qi_idx,' vs ',qi_test,'=qi_test'
					end if
				end do
			end do
		end do
		!
	end subroutine


	integer function get_rel_kpt(qix, qiy, qiz, kpt)
		integer,	intent(in)	::	qix, qiy, qiz
		real(dp),	intent(out)	::	kpt(3)
		!
		kpt(1)	=	(	 2.0_dp*real(qix,dp)	- real(mp_grid(1),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(1),dp) )
		kpt(2)	=	(	 2.0_dp*real(qiy,dp)	- real(mp_grid(2),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(2),dp) )
		kpt(3)	=	(	 2.0_dp*real(qiz,dp)	- real(mp_grid(3),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(3),dp) )
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