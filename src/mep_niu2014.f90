module mep_niu2014
	use constants,		only:		dp, i_dp 
	use input_paras,	only:		kubo_tol,	valence_bands
	use matrix_math,	only:		my_Levi_Civita,			& 
									convert_tens_to_vect

	implicit none

	private
	public					::		mep_niu2014_full
				
	save

contains

		!
!MEP RESPONSES
	pure function mep_niu2014_full(velo, en) result(F_ij)
		complex(dp),	intent(in)		::	velo(:,:,:)
		real(dp),		intent(in)		::	en(:)
		real(dp)						::	F_ij(3,3)
		real(dp)						::	dE_0n
		complex(dp)						::	w_n0(3)
		integer							::	n, n0, i
		!
		F_ij 	=	0.0_dp
		!
		do n0 =1 , valence_bands
			do n = 1, size(en)
				!
				if( n/= n0) then
					!
					w_n0(:)	=	get_w_n0(n, n0, velo, en)
					dE_0n	=	(	en(n0) - en(n)	)**2
					!
					do i = 1, 3
						F_ij(i,:)	=	F_ij(i,:)	+ real(		velo(i,n0,n) * w_n0(:) / dE_0n		, dp)
					end do
				end if
				!
			end do
		end do
		!
		return
	end function



	pure function get_w_n0(n, n0, velo, en) result(w_n0)
		integer,			intent(in)		::	n, n0
		complex(dp),		intent(in)		::	velo(:,:,:)
		real(dp),			intent(in)		::	en(:)
		complex(dp)							::	w_n0(3)
		complex(dp)							::	pre_fact
		real(dp)							::	dE_m0
		integer								::	m, j, k, l
		!
		w_n0	=	cmplx(	0.0_dp	,	0.0_dp,			dp)
		!
		!
		do j = 1, 3
			do k = 1,3 
				do l = 1, 3
					pre_fact	= 	- i_dp * real(my_Levi_Civita(j,k,l),dp)
					!
					!
					do m = 1, size(en)
						if( m /= n0 ) then
							dE_m0	=	en(m)	- en(n0)
							!
							w_n0(j)	=	w_n0(j)	+ 	pre_fact	*	velo(k,n,m)	* velo(l,m,n0)	/ dE_m0
							!
							if( m == n ) then
								w_n0(j)	=	w_n0(j)		+  pre_fact		*	velo(k,n0,n0) * velo(l,n,n0)	/ dE_m0
							end if
						end if
					end do
					!
					!
				end do
			end do
		end do
		return
	end function




end module mep_niu2014


