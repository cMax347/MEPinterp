module photo

	use constants,		only:		dp, i_dp, pi_dp
	use statistics,		only:		fd_stat

	implicit none

	private
	public					::		photo_2nd_cond


contains


	pure function photo_2nd_cond(en_k, V_ka, hw_lst, phi_laser, fd_distrib, i_eta_smr) result(phot_cond)
		!
		!	implementation of Yang et al., PRB 97 241118(R) (2018)    EQ(1)
		!	"PHOTOGALVANIC EFFECT IN WEYL SEMIMETALS FROM FIRST PRINCIPLES"
		!
		!		EQ.(1)		
		!
		real(dp),		intent(in)			::	hw_lst(:), fd_distrib(:,:), en_k(:)
		complex(dp),	intent(in)			::	V_ka(:,:,:), phi_laser, i_eta_smr
		real(dp),		allocatable			::	phot_cond(:,:,:,:,:), tmp(:,:,:,:), df_ln(:)
		complex(dp)							::	dE_nm, dE_nl,  dE_nl_hw, v_nl_lm_nm(3,3,3)
		integer								::	m, n, l, c, b, hw, omega, ef_idx, n_ef, n_hw, n_wf
		!
		n_ef	=	size(	fd_distrib	,	1)
		n_hw 	=	size(	hw_lst		,	1)
		n_wf	=	size(	en_k		,	1)


		allocate(	df_ln(							n_ef	))
		allocate(	tmp(			3,3,3,	n_hw			))
		allocate(	phot_cond(		3,3,3,	n_hw,	n_ef	))
		phot_cond	=	0.0_dp
		!
		!	LOOP STATES
		do 	m = 1, n_wf
			do n = 1, n_wf
				dE_nm	=	cmplx(	en_k(n) - en_k(m)	,0.0_dp,dp) 	-	 i_eta_smr
				do l = 1, n_wf
					!
					df_ln	=	fd_distrib(:,l)	-	fd_distrib(:,n)
					dE_nl	=		en_k(n) 	- 		en_k(l)
					!
					!	LOOP DIRECTIONS
					do c = 1, 3
						do b = 1, 3
							v_nl_lm_nm(:,b,c)	=	V_ka(:,n,l) * (V_ka(b,l,m) * V_ka(c,m,n))
						end do
					end do
					!
					!	LOOP FREQUENCIES
					tmp =	0.0_dp
					do hw = 1, n_hw
						do omega = -1, 1, 2
							dE_nl_hw	=	dE_nl +	omega*hw_lst(hw) -	 i_eta_smr		
							!
							tmp(:,:,:,hw)	=	tmp(:,:,:,hw)	+	real(		phi_laser	*	v_nl_lm_nm(:,:,:)				&		
																					/	( dE_nm * dE_nl_hw * hw_lst(hw)**2	) 	,dp)			
						end do
					end do
					!
					!	LOOP FERMI LEVEL
					do ef_idx = 1, n_ef
						phot_cond(:,:,:,:,ef_idx)		=	tmp(:,:,:,:)	* df_ln(ef_idx)
					end do
					!
				end do
			end do
		end do
		!
		!
		return
	end function

end module photo