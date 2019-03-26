module photo

	use constants,		only:		dp, i_dp, pi_dp
	use statistics,		only:		fd_stat
	use omp_lib

	implicit none

	private
	public					::		photo_2nd_cond


contains


	function photo_2nd_cond(en_k, V_ka, hw_lst, phi_laser, fd_distrib, i_eta_smr) result(phot_cond)
		!
		!	implementation of Yang et al., PRB 97 241118(R) (2018)    EQ(1)
		!	"PHOTOGALVANIC EFFECT IN WEYL SEMIMETALS FROM FIRST PRINCIPLES"
		!
		!		EQ.(1)		
		!
		real(dp),		intent(in)			::	hw_lst(:), fd_distrib(:,:), en_k(:)
		complex(dp),	intent(in)			::	V_ka(:,:,:), phi_laser, i_eta_smr
		real(dp),		allocatable			::	phot_cond(:,:,:,:,:), df_ln(:)
		complex(dp)							::	dE_nm, dE_nl,  dE_nl_hw, vvv_nl_lm_mn(3,3,3)
		integer								::	m, n, l, c, b, hw, omega, ef_idx, n_ef, n_hw, n_wf
		!
		n_ef	=	size(	fd_distrib	,	1)
		n_hw 	=	size(	hw_lst		,	1)
		n_wf	=	size(	en_k		,	1)
		!
		allocate(	df_ln(							n_ef	))
		allocate(	phot_cond(		3,3,3,	n_hw,	n_ef	))
		phot_cond	=	0.0_dp
		!
		!
		!$OMP  PARALLEL DO DEFAULT(NONE)  																&
		!$OMP PRIVATE(n, m, l, dE_nm, df_ln, dE_nl, c, b, vvv_nl_lm_mn,  ef_idx,hw, omega, dE_nl_hw)	&
		!$OMP SHARED(n_wf, n_hw, n_ef, en_k, fd_distrib, V_ka, hw_lst, phi_laser, i_eta_smr) 			&
		!$OMP REDUCTION(+: phot_cond)
		do 	n = 1, n_wf
			do m = 1, n_wf
				dE_nm	=	cmplx(	en_k(n) - en_k(m)	,0.0_dp,dp) 			-	 i_eta_smr
				do l = 1, n_wf
					if(l==n) cycle
					!
					df_ln(:)=	fd_distrib(:,l)	-	fd_distrib(:,n)
					dE_nl	=	cmplx(	en_k(n) - 	 en_k(l),	0.0_dp, dp)		-	 i_eta_smr
					!
					!	LOOP DIRECTIONS
					do c = 1, 3
						do b = 1, 3
							vvv_nl_lm_mn(:,b,c)	=	V_ka(:,n,l) * (V_ka(b,l,m) * V_ka(c,m,n))
						end do
					end do
					!
					!	LOOP FERMI LEVEL
					do ef_idx =1, n_ef
						if(abs(df_ln(ef_idx))>1e-6_dp) then
							!	LOOP FREQUENCIES
							do hw = 1, n_hw
								do omega = -1, 1, 2
									dE_nl_hw	=	dE_nl +	cmplx(	real(omega,dp)*hw_lst(hw),	0.0_dp,	dp) 		
									!
									phot_cond(:,:,:,hw,ef_idx)	=	phot_cond(:,:,:,hw,ef_idx)							&
										+	real(	phi_laser *	vvv_nl_lm_mn(:,:,:)	/ ( dE_nm * dE_nl_hw ) 		,dp)	&
										/ 	hw_lst(hw)**2		
								end do
							end do
							phot_cond(:,:,:,:,ef_idx)		=	phot_cond(:,:,:,:,ef_idx)	*	df_ln(ef_idx)
						end if
					end do 
					!
					!
				end do
			end do
		end do
		!$OMP END PARALLEL DO
		!
		!
		return
	end function

end module photo