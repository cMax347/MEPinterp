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
		real(dp)							::	dE_nm, dE_nl, re_dE_smr,re_dE,im_dE, im_dE_smr
		complex(dp)							::	dE_nm_smr, dE_nl_smr, dE_nl_hw, vvv_nl_lm_mn(3,3,3)
		integer								::	m, n, l, a, b, hw,smr, omega, ef_idx, n_ef, n_hw, n_wf, n_smr
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
		!$OMP  PARALLEL DO DEFAULT(PRIVATE) COLLAPSE(2) 																&
		!!!!$OMP PRIVATE(l, dE_nm,dE_nm_smr, df_ln,re_dE_denom,im_dE_denom,cmplx_dE_denom, dE_nl,dE_nl_smr,dE_nl_hw, a, b, vvv_nl_lm_mn, tmp, hw, omega, ef_idx, smr)	&
		!$OMP SHARED(n_wf, n_hw,n_smr, n_ef, en_k, fd_distrib, V_ka, hw_lst, eta_smr_lst) 							&
		!$OMP REDUCTION(+: phot_cond)
		do 	n = 1, n_wf
			do m = 1, n_wf
				!if(m==n) cycle
				!dE_nm	=	en_k(n) - en_k(m)
				do l = 1, n_wf
					if(l==n) cycle
					!	ZERO INIT
					tmp 			=	cmplx(0.0_dp,0.0_dp,dp)
					vvv_nl_lm_mn 	= 	cmplx(0.0_dp,0.0_dp,dp)
					!
					!	ENERGY DENOM
					re_dE	=	  en_k(n)**2 		- en_k(n)*en_k(l) 	- en_k(m)*en_k(n)	+ en_k(m)*en_k(l) 
					im_dE	= 	-2.0_dp * en_k(n) 	+ en_k(m) 			+ en_k(l)
					!
					!	LOOP DIRECTIONS
					do a = 1, 3
						do b = 1, 3
							vvv_nl_lm_mn(:,a,b)	=	V_ka(a,n,l) * V_ka(b,l,m) * V_ka(:,m,n)
						end do
					end do
					!
					!	LOOP SMEARING
					do smr = 1, n_smr
						!	LOOP FREQUENCIES
						do hw = 1, n_hw
							do omega = -1, 1, 2
								!dE_nl_hw	=	dE_nl_smr	+ real(omega,dp)	* hw_lst(hw)
								!!if(smr==1) &
								!!write(*,'(a,i1,a,i1,a,i1,a,f6.2,a,f6.3,a,20(e12.4))')	'dE-DENOM_',n,',',l,',',m,'_smr(omega=',real(omega,dp),&
								!!												"; hw=",hw_lst(hw),")=",1.0_dp	/	(dE_nm_smr*dE_nl_hw)
								!!
								!tmp(:,:,:,hw,smr)	=	tmp(:,:,:,hw,smr)	+	vvv_nl_lm_mn(:,:,:)				&		
								!												/	(dE_nm_smr*dE_nl_hw)
								!					!	
								!					!	/	( dE_nm_smr * (dE_nl_smr +	cmplx(omega,0.0_dp,dp)*hw_lst(hw)) * hw_lst(hw)**2	) 
								!
								!
								re_dE_smr			=	re_dE   								& 
														+ real(omega,dp)*hw_lst(hw)*en_k(n)		& 
														- real(omega,dp)*hw_lst(hw)*en_k(m) 	&
														- eta_smr_lst(smr)**2				
								!
								im_dE_smr			=	eta_smr_lst(smr) * ( im_dE - real(omega,dp)*hw_lst(hw)	)
								!
								tmp(:,:,:,hw,smr)	=	tmp(:,:,:,hw,smr)	+	vvv_nl_lm_mn(:,:,:)				&		
																				/	cmplx(re_dE_smr, im_dE_smr,dp)
								!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~				
							end do
							tmp(:,:,:,hw,smr)		=	tmp(:,:,:,hw,smr)	/	(hw_lst(hw)*hw_lst(hw))
						end do
						!
					end do
					!
					!
					!	LOOP FERMI LEVEL
					df_ln(:)	=	fd_distrib(:,l)	-	fd_distrib(:,n)
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
		return
	end function

end module photo