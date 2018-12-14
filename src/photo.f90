module photo

	use constants,		only:		dp, i_dp, pi_dp
	use input_paras,	only:		kubo_tol
	use statistics,		only:		fd_stat

	implicit none

	private
	public					::		photo_2nd_cond


contains


	pure function photo_2nd_cond(hw, phi_laser, e_fermi, T_kelvin, i_eta_smr, en_k, V_ka) result(phot_cond)
		!
		!	implementation of Yang et al., PRB 97 241118(R) (2018) 
		!	"PHOTOGALVANIC EFFECT IN WEYL SEMIMETALS FROM FIRST PRINCIPLES"
		!
		!		EQ.(1)		
		!
		real(dp),		intent(in)			::	hw, e_fermi, T_kelvin, en_k(:)
		complex(dp),	intent(in)			::	V_ka(:,:,:), phi_laser, i_eta_smr
		real(dp)							::	phot_cond(3,3)
		complex(dp)							::	cmplx_cond(3,3), dE_nm, dE_nl, df_ln
		real(dp)							:: 	pre_fact, w
		integer								::	m, n, l, omega, b, c
		!
		phot_cond	=	0.0_dp
		!
		if( abs(hw)	> kubo_tol ) then
			!	IN ATOMIC UNITS:
			w			=	hw		
			pre_fact	=	1.0_dp / 	(  w**2)
			!
			!
			cmplx_cond	=	cmplx(0.0_dp,0.0_dp,dp)
			do 	m = 1, size(en_k)
				do n = 1, size(en_k)
					dE_nm	=	cmplx(	en_k(n) - en_k(m)	,0.0_dp,dp) 	-	 i_eta_smr
					!
					do l = 1, size(en_k)
						df_ln	=	fd_stat(en_k(l),e_fermi,T_kelvin)	-	fd_stat(en_k(n),e_fermi, T_kelvin)	
						!
						do omega = -1, 1, 2
							dE_nl	=	cmplx(	en_k(n) - en_k(l) +	omega*hw	,0.0_dp,dp) 	-	 i_eta_smr
							!
							do b = 1, 3
								do c = 1, 3
									cmplx_cond(:,b)	=	cmplx_cond(:,b)		+	df_ln *		( V_ka(:,n,l) * V_ka(b,l,m) * V_ka(c,m,n) )	/	( dE_nm * dE_nl	) 
								end do
							end do
						end do
					end do
				end do
			end do
			!
			cmplx_cond	=	phi_laser	*		cmplx_cond
			phot_cond	=	dreal(				cmplx_cond	)	/ 	hw**2
		end if
		!
		return
	end function

end module photo