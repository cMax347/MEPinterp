module photo

	use constants,		only:		dp, i_dp, pi_dp
	use statistics,		only:		fd_stat
	use omp_lib

	implicit none

	private
	public					::		photo_2nd_cond


contains


	function photo_2nd_cond(en_k, V_ka, hw_lst,eta_smr_lst, fd_distrib) result(phot_cond)
		!
		!	implementation of Yang et al., PRB 97 241118(R) (2018)    EQ(1)
		!	"PHOTOGALVANIC EFFECT IN WEYL SEMIMETALS FROM FIRST PRINCIPLES"
		!
		!		EQ.(1)		
		!
		real(dp),		intent(in)			::	hw_lst(:), eta_smr_lst(:), fd_distrib(:,:), en_k(:)
		complex(dp),	intent(in)			::	V_ka(:,:,:)
		real(dp),		allocatable			::	df_ln(:)
		complex(dp),	allocatable			::	phot_cond(:,:,:,:,:,:), tmp(:,:,:,:,:)
		real(dp)							::	dE_nm, dE_nl
		complex(dp)							::	dE_nm_smr, dE_nl_smr,  dE_nl_hw, vvv_nl_lm_mn(3,3,3), i_eta_smr
		integer								::	m, n, l, c, b, hw,smr, omega, ef_idx, n_ef, n_hw, n_wf, n_smr
		!
		n_ef	=	size(	fd_distrib	,	1)
		n_hw 	=	size(	hw_lst		,	1)
		n_wf	=	size(	en_k		,	1)
		n_smr	=	size(	eta_smr_lst ,	1)
		!
		allocate(	df_ln(									n_ef	))
		allocate(	tmp(			3,3,3,	n_hw,	n_smr			))
		allocate(	phot_cond(		3,3,3,	n_hw,	n_smr, 	n_ef	))
		phot_cond	=	0.0_dp
		!
		!	
		!$OMP  PARALLEL DO DEFAULT(NONE) COLLAPSE(2) 																&
		!$OMP PRIVATE(l, dE_nm,dE_nm_smr, df_ln, dE_nl,dE_nl_smr, c, b, vvv_nl_lm_mn, tmp, hw, omega, dE_nl_hw, ef_idx, smr, i_eta_smr)	&
		!$OMP SHARED(n_wf, n_hw,n_smr, n_ef, en_k, fd_distrib, V_ka, hw_lst, eta_smr_lst) 							&
		!$OMP REDUCTION(+: phot_cond)
		do 	n = 1, n_wf
			do m = 1, n_wf
				dE_nm	=	en_k(n) - en_k(m)
				do l = 1, n_wf
					if(l==n) cycle
					dE_nl 	=	en_k(n) 		- 	en_k(l)
					df_ln(:)=	fd_distrib(:,l)	-	fd_distrib(:,n)
					!
					!	LOOP DIRECTIONS
					do c = 1, 3
						do b = 1, 3
							vvv_nl_lm_mn(:,b,c)	=	V_ka(:,n,l) * (V_ka(b,l,m) * V_ka(c,m,n))
						end do
					end do
					!
					!	LOOP SMEARING
					do smr = 1, n_smr
						dE_nm_smr	=	cmplx(	dE_nm,	-eta_smr_lst(smr), dp) 	
						dE_nl_smr	=	cmplx(	dE_nl,	-eta_smr_lst(smr), dp)		
						!
						!	LOOP FREQUENCIES
						tmp =	0.0_dp
						do hw = 1, n_hw
							do omega = -1, 1, 2
								dE_nl_hw	=	dE_nl_smr +	real(omega,dp)*hw_lst(hw)	
								!
								tmp(:,:,:,hw,smr)	=	tmp(:,:,:,hw,smr)	+	vvv_nl_lm_mn(:,:,:)				&		
														/	( dE_nm_smr * dE_nl_hw * hw_lst(hw)**2	) 		
							end do
						end do
						!
					end do
					!
					!	LOOP FERMI LEVEL
					do ef_idx = 1, n_ef
						phot_cond(:,:,:,:,:,ef_idx)		=	phot_cond(:,:,:,:,:,ef_idx)	+	tmp(:,:,:,:,:)	* df_ln(ef_idx)
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