module photo

	use constants,		only:		dp, i_dp, pi_dp
	use input_paras,	only:		kubo_tol
	use statistics,		only:		fd_stat

	implicit none

	private
	public					::		photo_2nd_cond


contains


	pure function photo_2nd_cond(hw_lst, phi_laser, e_fermi, T_kelvin, i_eta_smr, en_k, V_ka) result(phot_cond)
		!
		!	implementation of Yang et al., PRB 97 241118(R) (2018)    EQ(1)
		!	"PHOTOGALVANIC EFFECT IN WEYL SEMIMETALS FROM FIRST PRINCIPLES"
		!
		!		EQ.(1)		
		!
		real(dp),		intent(in)			::	hw_lst(:), e_fermi, T_kelvin, en_k(:)
		complex(dp),	intent(in)			::	V_ka(:,:,:), phi_laser, i_eta_smr
		real(dp),		allocatable			::	phot_cond(:,:,:,:)
		complex(dp)							::	dE_nm, dE_nl, df_ln, dE_nl_hw, v_nl_lm_nm(3,3,3)
		integer								::	m, n, l, c, b, hw, omega
		!
		allocate(phot_cond(	3,3,3,size(hw_lst,1)))
		phot_cond	=	0.0_dp
		!
		!	LOOP STATES
		do 	m = 1, size(en_k)
			do n = 1, size(en_k)
				dE_nm	=	cmplx(	en_k(n) - en_k(m)	,0.0_dp,dp) 	-	 i_eta_smr
				do l = 1, size(en_k)
					dE_nl	=	en_k(n) - en_k(l)
					df_ln	=	fd_stat(en_k(l),e_fermi,T_kelvin)	-	fd_stat(en_k(n),e_fermi, T_kelvin)	
					!
					!	LOOP DIRECTIONS
					do c = 1, 3
						do b = 1, 3
							v_nl_lm_nm(:,b,c)	=	V_ka(:,n,l) * (V_ka(b,l,m) * V_ka(c,m,n))
						end do
					end do
					v_nl_lm_nm			=	v_nl_lm_nm	* df_ln
					!
					!
					!	LOOP FREQUENCIES
					do hw = 1, size(hw_lst,1)
						do omega = -1, 1, 2
							dE_nl_hw	=	dE_nl +	omega*hw_lst(hw) -	 i_eta_smr		
							!
							phot_cond(:,:,:,hw)	=	phot_cond(:,:,:,hw)	+	real(		v_nl_lm_nm(:,:,:)							&		
																					/	( dE_nm * dE_nl_hw * hw_lst(hw)**2	) 	,dp)			
						end do
					end do
					!
					!
				end do
			end do
		end do
		!
		!
		return
	end function

end module photo