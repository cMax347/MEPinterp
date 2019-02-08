module kubo
	!
	!	calculates some of the Kubo formulas available in postw90
	!	see wann guide chapter 12 for details
	!		http://www.wannier.org/doc/user_guide.pdf
	!
	!
	use omp_lib
	use input_paras,	only:	kubo_tol
	use constants,		only:	dp, i_dp, pi_dp, aUtoEv
	use statistics,		only:	fd_stat


	implicit none


	private
	public			::			kubo_ahc_tens,		&
								velo_ahc_tens,		&
								kubo_opt_tens,		&
								kubo_ohc_tens	



contains





!
!	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				ANOMALOUS HALL (hw -> 0 LIMIT)!	------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!public

	function kubo_ahc_tens(en_k, Om_kab, fd_distrib) result( o_ahc)
		!	see wann guide chapter 12
		!		eq. (12.16)
		real(dp),						intent(in)		::	en_k(:), fd_distrib(:,:)
		complex(dp), 	allocatable,	intent(in)		::	Om_kab(:,:,:,:)
		real(dp)										::	curv_nn(3,3)
		real(dp),		allocatable						::	o_ahc(:,:,:)
		integer											::	n, ef_idx, n_ef, n_wf
		!
		n_wf	=	size(en_k,1)
		n_ef	=	size(fd_distrib,1)
		allocate(		o_ahc(	3,3		,	n_ef	))
		o_ahc	=	0.0_dp
		!
		!
		if(allocated(Om_kab)) then
			!$OMP PARALLEL  DO 									&
			!$OMP DEFAULT(	none							)	&	
			!$OMP PRIVATE( 	curv_nn, ef_idx					)	&
			!$OMP SHARED(	n_wf, n_ef, Om_kab, fd_distrib	)	&
			!----
			!$OMP 	REDUCTION( + : o_ahc	)
			do n = 1, n_wf
				curv_nn(:,:)	=	real(Om_kab(:,:,n,n),dp)
				!
				do ef_idx =1 ,n_ef
					o_ahc(:,:,ef_idx)	=	o_ahc(:,:,ef_idx)	-	curv_nn(:,:)	*	fd_distrib(ef_idx,n)
				end do
			end do
			!$OMP END PARALLEL DO
		end if
		!
		!
		return
	end function

!
!
!	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				OPTICAL CONDUTCTIVITY (hw \non_eq 0)
!	------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	function velo_ahc_tens(en_k, v_kab, hw_lst, fd_distrib, i_eta_smr) result( o_ahc)
		!
		!	use Wanxiangs expression to calculate the ohc
		!
		!		VIA VELOCITIES
		real(dp),		intent(in)		::	en_k(:), hw_lst(:), fd_distrib(:,:)
		complex(dp), 	intent(in)		::	v_kab(:,:,:), i_eta_smr
		complex(dp)	,	allocatable		::	o_ahc(:,:,:,:), vv_dE_hw(:,:,:)
		real(dp),		allocatable		::	d_ef_lst(:)
		real(dp)						::	dE_sqr, v_nm_mn(3,3)
		integer							::	n, l, j, hw_idx, n_ef, n_hw, n_wf, ef_idx
		!
		n_hw	=	size(	hw_lst		,	1)
		n_ef	=	size(	fd_distrib	,	1)
		n_wf	=	size(	en_k		,	1)
		!
		allocate(		d_ef_lst(						n_ef	)	)
		allocate(		vv_dE_hw(	3,3		,	n_hw			)	)
		allocate(	o_ahc(			3,3		,	n_hw,	n_ef	)	)
		o_ahc	=	0.0_dp
		!
		!$OMP  PARALLEL DO DEFAULT(NONE)  											&
		!$OMP PRIVATE(l, d_ef_lst, dE_sqr, j, v_nm_mn, hw_idx, vv_dE_hw, ef_idx)	&
		!$OMP SHARED(fd_distrib, en_k, v_kab, hw_lst, n_wf, n_hw, n_ef, i_eta_smr) 	&
		!$OMP REDUCTION(+: o_ahc)
		do n = 1, n_wf
			do l = 1, n_wf
				if(l/=n) then
					!
					!	dE	BANDS
					d_ef_lst	=	fd_distrib(:,n)		-	fd_distrib(:,l)
					dE_sqr		= 	(	en_k(n) 		- 		en_k(l)	  )**2 	
					!
					!	loop real space direction
					do j = 1, 3
						v_nm_mn(:,j)			=	aimag(	v_kab(:,n,l) * v_kab(j,l,n)		) 
					end do
					!
					!	LOOP HW
					do hw_idx =1, n_hw
						vv_dE_hw(:,:,hw_idx)	=	 v_nm_mn(:,:)	/	(	dE_sqr	- (	hw_lst(hw_idx) + i_eta_smr	)**2	)	
					end do
					!
					!	LOOP FERMI LEVEL
					do ef_idx = 1, n_ef
						o_ahc(:,:,:,ef_idx)		=	o_ahc(:,:,:,ef_idx)	+	d_ef_lst(ef_idx)	*	vv_dE_hw(:,:,:)
					end do
				end if
				!
				!
			end do
		end do
		!$OMP END PARALLEL DO
		!
		!
		return
	end function




	function kubo_ohc_tens(en_k, v_kab, hw_lst, fd_distrib, i_eta_smr)	result(z_ohc)
		!
		!		calculates wann guide (12.5) explicitly
		!	
		!	via VELOCITIES
		real(dp),		intent(in)		::	en_k(:), hw_lst(:), fd_distrib(:,:)
		complex(dp),	intent(in)		::	v_kab(:,:,:), i_eta_smr
		complex(dp), 	allocatable		::	z_ohc(:,:,:,:), vv_dE_hw(:,:,:) 
		real(dp), 		allocatable		::	df_mn(:) 
		complex(dp)						::	v_nm_mn(3,3)
		real(dp)						::	dE_mn
		integer							::	n, m, j, n_ef, n_hw, n_wf, ef_idx, hw_idx
		!
		n_ef	=	size(	fd_distrib	,	1	)
		n_hw	=	size(	hw_lst		,	1	)
		n_wf	=	size(	en_k		,	1	)
		!
		allocate(		df_mn(						n_ef	))
		allocate(		vv_dE_hw(	3,3,	n_hw			))
		allocate(		z_ohc(		3,3,	n_hw,	n_ef	))
		!
		z_ohc	=	cmplx(0.0_dp, 0.0_dp, dp)
		!
		!$OMP  PARALLEL DO DEFAULT(NONE)  													&
		!$OMP PRIVATE(n, df_mn, dE_mn, j, v_nm_mn, vv_dE_hw, hw_idx, ef_idx)				&
		!$OMP SHARED(fd_distrib, en_k, v_kab, hw_lst, n_wf, n_hw, n_ef, kubo_tol, i_eta_smr)&
		!$OMP REDUCTION(+: z_ohc)
		do m = 1 , n_wf
			do n = 1, n_wf
				dE_mn 	=	en_k(m)	-	en_k(n)
				!
				!
				if(	abs(dE_mn) > kubo_tol ) then 
					df_mn		=	fd_distrib(:,m)	-	fd_distrib(:,n)
					!
					!	loop real space direction
					do j = 1, 3
						v_nm_mn(:,j)	=	v_kab(:,n,m) * v_kab(j,m,n)
					end do
					!
					!	LOOP HW
					vv_dE_hw	=	cmplx(0.0_dp,0.0_dp,dp)
					do hw_idx = 1, n_hw
						vv_dE_hw(:,:,hw_idx)	=	v_nm_mn(:,:)	/	(	dE_mn * (dE_mn - hw_lst(hw_idx) - i_eta_smr))
					end do
					!
					!	LOOP FERMI LEVEL
					do ef_idx = 1, n_ef
						z_ohc(:,:,:,ef_idx)		=	z_ohc(:,:,:,ef_idx)	+ 	i_dp * 	vv_dE_hw(:,:,:) * df_mn(ef_idx)
					end do
				end if
			end do
		end do
		!$OMP END PARALLEL DO
		!
		return
	end function





	subroutine kubo_opt_tens(en_k, A_ka, hw_lst, fd_distrib, i_eta_smr, opt_symm, opt_asymm)	
		!	see wann guide chapter 12
		!		eq. (12.14) & (12.15)
		!
		!	via CONNECTION
		real(dp),					intent(in)		::	en_k(:), hw_lst(:), fd_distrib(:,:) 
		complex(dp),				intent(in)		::	i_eta_smr
		complex(dp), allocatable,	intent(in)		::	A_ka(:,:,:)
		complex(dp), allocatable,	intent(out)		::	opt_symm(:,:,:,:), 	opt_asymm(:,:,:,:)
		complex(dp), allocatable					::	opt_herm(:,:,:,:),	opt_aherm(:,:,:,:)
		integer										::	n_hw, n_ef
		!
		n_hw	=	size(	hw_lst		,	1)
		n_ef	=	size(	fd_distrib	,	1)
		!
		allocate(	opt_symm(	3,3,	n_hw, n_ef		))
		allocate(	opt_asymm(	3,3,	n_hw, n_ef		))
		!
		opt_symm	= 	0.0_dp
		opt_asymm	=	0.0_dp
		!
		if(allocated(A_ka)) then
			allocate(	opt_aherm(	3,3,	n_hw, n_ef		))
			allocate(	opt_herm(	3,3,	n_hw, n_ef		))
			!
			opt_herm	=	get_hermitian(hw_lst, fd_distrib, i_eta_smr, en_k, A_ka)
			opt_aherm	=	get_anti_herm(hw_lst, fd_distrib, i_eta_smr, en_k, A_ka)
			!
			opt_symm	=	real(	opt_herm	,	dp) 		+	i_dp * aimag(	opt_aherm	)		!	(12.14)		
			opt_asymm	=	real(	opt_aherm	,	dp)			+	i_dp * aimag(	opt_herm	)		!	(12.15)
		end if
		!
		return
	end subroutine









!private:
	real(dp) pure function delta_broad(eps, smr)
		!see wannier guide eq. (12.9)
		! 'broadened delta-function'
		real(dp),		intent(in)	::	eps		!	energy
		complex(dp),	intent(in)	::	smr		!	smearing factor
		!
		if(abs(eps-smr) > 1e-6_dp)	then
			delta_broad		=	aimag(	1.0_dp/ ( eps - smr )			)			/	pi_dp
		else
			delta_broad		= 	aimag(1.0_dp / (1e-6_dp*i_dp)	)	/	pi_dp		!
		end if
		return
	end function



	function get_hermitian(hw_lst, fd_distrib, i_eta_smr, en_k, A_ka)	result(opt_herm)
		!	wann guide eq. (12.10)
		real(dp),		intent(in)		::	hw_lst(:), fd_distrib(:,:), en_k(:)
		complex(dp),	intent(in)		::	i_eta_smr
		complex(dp),	intent(in)		::	A_ka(:,:,:)
		complex(dp),	allocatable		::	opt_herm(:,:,:,:), a_dE_hw(:,:,:)
		real(dp),		allocatable		::	df_mn(:)
		complex(dp)						::	a_nm_mn(3,3)
		integer							::	n, m, b, 	&
											hw_idx, ef_idx, n_hw, n_ef, n_wf 
		real(dp)						::	dE_mn
		!
		n_hw	=	size(hw_lst		,	1)
		n_ef	=	size(fd_distrib	,	1)
		n_wf	=	size(en_k		,	1)
		!
		allocate(	df_mn(						n_ef	))
		allocate(	a_dE_hw(	3,3,	n_hw			))
		allocate(	opt_herm(	3,3,	n_hw,	n_ef	))

		opt_herm	=	0.0_dp
		!
		do m = 1, n_wf
			do n = 1, n_wf
				if(n/=m) then
					!
					!	dE	BANDS				
					df_mn(:)	=	fd_distrib(:,m)	-	fd_distrib(:,n)
					dE_mn		=		en_k(m)		-	en_k(n) 
					!
					!	
					do b = 1, 3
						a_nm_mn(:,b)	=	A_ka(:,n,m) * A_ka(b,m,n)
					end do
					a_nm_mn				=	dE_mn		* a_nm_mn
					!
					!	LOOP HW
					a_dE_hw	=	cmplx(0.0_dp,0.0_dp,dp)
					do hw_idx = 1, n_hw
						a_dE_hw(:,:,hw_idx)	=	a_nm_mn(:,:)			&
												* delta_broad( dE_mn - hw_lst(hw_idx)	, i_eta_smr)
					end do
					!
					!	LOOP FERMI LEVEL
					do ef_idx = 1, n_ef
						opt_herm(:,:,:,ef_idx)	=	opt_herm(:,:,:,ef_idx)	+	a_dE_hw(:,:,:)	* df_mn(ef_idx)
					end do
					!
				end if
			end do
		end do
		!
		! 
		opt_herm	=	cmplx(-	pi_dp,0.0_dp, dp)	* opt_herm
		!
		return
	end function



	function get_anti_herm(hw_lst, fd_distrib, i_eta_smr, en_k, A_ka)	result(opt_aherm)
		!	wann guide eq. (12.11)
		real(dp),		intent(in)		::	hw_lst(:), fd_distrib(:,:), en_k(:)
		complex(dp),	intent(in)		::	i_eta_smr	
		complex(dp),	intent(in)		::	A_ka(:,:,:)
		complex(dp),	allocatable		::	a_dE_hw(:,:,:) ,opt_aherm(:,:,:,:)
		real(dp),		allocatable		::	df_mn(:)
		complex(dp)						::	a_nm_mn(3,3)
		integer							::	n, m, b, 		&
											hw_idx, n_hw,	&
											ef_idx, n_ef,	&
											n_wf 
		real(dp)						::	dE_mn 
		!
		n_hw	=	size(	hw_lst		,	1	)
		n_ef	=	size(	fd_distrib	,	1	)
		n_wf	=	size(	en_k		,	1	)
		!
		allocate(	df_mn(						n_ef	))
		allocate(	a_dE_hw(	3,3,	n_hw			))
		allocate(	opt_aherm(	3,3,	n_hw,	n_ef	))
		opt_aherm	=	cmplx(0.0, 0.0, dp)
		!
		do m = 1, n_wf
			do n = 1, n_wf
				if(n/=m) then
					!	dE	BANDS
					df_mn(:)	=	fd_distrib(:,m)		-	fd_distrib(:,n)
					dE_mn		=		en_k(m) 		- 		en_k(n)	
					!
					do b = 1, 3
						a_nm_mn(:,b)	=	A_ka(:,n,m) * A_ka(b,m,n)	
					end do
					a_nm_mn		=	a_nm_mn  * dE_mn
					!
					!	LOOP HW LASER
					a_dE_hw		=	cmplx(0.0_dp,0.0_dp,dp)
					do hw_idx = 1, n_hw
						a_dE_hw(:,:,hw_idx)	=	a_nm_mn(:,:)							&
												/	real(dE_mn - hw_lst(hw_idx) - i_eta_smr , dp)
					end do
					!
					!	LOOP FERMI LEVEL
					do ef_idx = 1, n_ef
						opt_aherm(:,:,:,ef_idx)	=	opt_aherm(:,:,:,ef_idx)	+ a_dE_hw(:,:,:)	* df_mn(ef_idx)
					end do
					!
				end if
			end do
		end do
		!
		opt_aherm	=	i_dp	*	opt_aherm
		!
		return
	end function


















































end module kubo