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
	public			::			kubo_DC_ahc,			&
								kubo_AC_ahc,			&
								kubo_opt_tens,			&
								kubo_ohc_tens	



contains





!
!	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				ANOMALOUS HALL (hw -> 0 LIMIT)!	------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!public

	function kubo_DC_ahc(en_k, V_ka, fd_distrib) result( o_ahc)
		!	see wann guide chapter 12
		!		eq. (12.16)
		real(dp),						intent(in)		::	en_k(:), fd_distrib(:,:)
		complex(dp), 					intent(in)		::	V_ka(:,:,:)
		real(dp)										::	curv_nn(3,3)
		real(dp),		allocatable						::	o_ahc(:,:,:)
		integer											::	n, np, j, ef_idx, n_ef, n_wf
		!
		n_wf	=	size(en_k,1)
		n_ef	=	size(fd_distrib,1)
		allocate(		o_ahc(	3,3		,	n_ef	))
		o_ahc	=	0.0_dp
		!
		!
		!$OMP PARALLEL  DO 								&
		!$OMP DEFAULT(none)								&	
		!$OMP PRIVATE(np, j, curv_nn, ef_idx)			&
		!$OMP SHARED(n_wf, n_ef,en_k, V_ka, fd_distrib)	&
		!----
		!$OMP 	REDUCTION( + : o_ahc	)
		do n = 1, n_wf
			!
			curv_nn	=	0.0_dp
			do np = 1, n_wf
				if(np==n)	cycle
				do j = 1, 3
					curv_nn(:,j)	=	curv_nn(:,j)	- 2.0_dp *	aimag(	V_ka(:,n,np) * V_ka(j,np,n)	)	/ (en_k(np) - en_k(n))**2
				end do
			end do
			!
			do ef_idx =1 ,n_ef
				o_ahc(:,:,ef_idx)	=	o_ahc(:,:,ef_idx)	-	curv_nn(:,:)	*	fd_distrib(ef_idx,n)
			end do
		end do
		!$OMP END PARALLEL DO
		!
		!
		return
	end function

!
!
!	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				OPTICAL CONDUTCTIVITY (hw \non_eq 0)
!	------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	function kubo_AC_ahc(en_k, V_ka, hw_lst, fd_distrib, i_eta_smr) result( o_ahc)
		!
		!	use Wanxiangs expression to calculate the ohc
		!
		!		VIA VELOCITIES
		real(dp),		intent(in)		::	en_k(:), hw_lst(:), fd_distrib(:,:)
		complex(dp), 	intent(in)		::	V_ka(:,:,:), i_eta_smr
		complex(dp)	,	allocatable		::	o_ahc(:,:,:,:), vv_dE_hw(:,:,:)
		real(dp),		allocatable		::	d_ef_lst(:)
		real(dp)						::	dE_sqr, v_nm_mn(3,3)
		integer							::	n, np, j, hw_idx, n_ef, n_hw, n_wf, ef_idx
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
		!$OMP PRIVATE(np, d_ef_lst, dE_sqr, j, v_nm_mn, hw_idx, vv_dE_hw, ef_idx)	&
		!$OMP SHARED(fd_distrib, en_k, V_ka, hw_lst, n_wf, n_hw, n_ef, i_eta_smr) 	&
		!$OMP REDUCTION(+: o_ahc)
		do n = 1, n_wf
			do np = 1, n_wf
				if(np==n) cycle
				!
				!	dE	BANDS
				d_ef_lst	=		fd_distrib(:,n)		-	fd_distrib(:,np)
				dE_sqr		= 	(	en_k(	 n) 		- 		en_k(	 np)	 )**2 	
				!
				!	loop real space direction
				do j = 1, 3
					v_nm_mn(:,j)			=	aimag(	V_ka(:,n,np) * V_ka(j,np,n)		) 
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
				!
			end do
		end do
		!$OMP END PARALLEL DO
		!
		!
		return
	end function



!	^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				OPTICAL CONDUTCTIVITY (postw90 version, DEPRECATED)
!	------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	function kubo_ohc_tens(en_k, V_ka, hw_lst, fd_distrib, i_eta_smr)	result(z_ohc)
		!
		!		calculates wann guide (12.5) explicitly
		!	
		!	via VELOCITIES
		real(dp),		intent(in)		::	en_k(:), hw_lst(:), fd_distrib(:,:)
		complex(dp),	intent(in)		::	V_ka(:,:,:), i_eta_smr
		complex(dp), 	allocatable		::	z_ohc(:,:,:,:), vv_dE_hw(:,:,:) 
		real(dp), 		allocatable		::	df_mn(:) 
		complex(dp)						::	v_nm_mn(3,3)
		real(dp)						::	dE_mn
		integer							::	n, m, j, n_ef, n_hw, n_wf, ef_idx, hw_idx
		!
		write(*,*)	"[kubo_ohc_tens]: WARNING THIS FUNCTION IS DEPRECATED (consider using kubo_AC_ahc function instead)"
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
		!$OMP SHARED(fd_distrib, en_k, V_ka, hw_lst, n_wf, n_hw, n_ef, kubo_tol, i_eta_smr)&
		!$OMP REDUCTION(+: z_ohc)
		do m = 1 , n_wf
			do n = 1, n_wf
				if(m==n)	cycle
				dE_mn 		=	en_k(m)			-	en_k(n)
				df_mn		=	fd_distrib(:,m)	-	fd_distrib(:,n)
				!
				!	loop real space direction
				do j = 1, 3
					v_nm_mn(:,j)	=	V_ka(:,n,m) * V_ka(j,m,n)
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
	 function delta_broad(dE_mn, hw_lst,  smr)	result(d_broad)
		!see wannier guide eq. (12.9)
		! 'broadened delta-function'
		real(dp),		intent(in)	::	dE_mn, hw_lst(:)		!	energy
		complex(dp),	intent(in)	::	smr				!	smearing factor
		complex(dp),	allocatable	::	d_broad(:)
		complex(dp)					::	dE_smr
		integer						::	n_hw, hw_idx
		!
		n_hw	=	size(hw_lst,1)
		allocate(	d_broad(n_hw))
		d_broad	=	cmplx(0.0_dp,0.0_dp,dp)
		!
		dE_smr	=	pi_dp * (dE_mn - smr )
		!
		do hw_idx = 1, n_hw
			d_broad(hw_idx)		=	aimag(	1.0_dp	/ 	(dE_smr	-	pi_dp * hw_lst(hw_idx))	)			
		end do
		!
		return
	end function



	function get_hermitian(hw_lst, fd_distrib, i_eta_smr, en_k, A_ka)	result(opt_herm)
		!	wann guide eq. (12.10)
		real(dp),		intent(in)		::	hw_lst(:), fd_distrib(:,:), en_k(:)
		complex(dp),	intent(in)		::	i_eta_smr
		complex(dp),	intent(in)		::	A_ka(:,:,:)
		complex(dp),	allocatable		::	opt_herm(:,:,:,:), a_dE_hw(:,:,:), d_broad(:)
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
		!
		opt_herm	=	0.0_dp
		!
		!
		!$OMP PARALLEL DO DEFAULT(none)												&
		!$OMP PRIVATE(df_mn, dE_mn, b, a_nm_mn, d_broad, a_dE_hw, hw_idx, ef_idx )	&
		!$OMP SHARED(n_wf, n_hw, n_ef,  hw_lst, fd_distrib, en_k, A_ka, i_eta_smr)	&
		!$OMP REDUCTION(+: opt_herm)
		do m = 1, n_wf
			do n = 1, n_wf
				if(n/=m) then
					!
					!	dE	BANDS				
					df_mn(:)				=	fd_distrib(:,m)	-	fd_distrib(:,n)
					dE_mn					=		en_k(m)		-	en_k(n) 
					! --
					!
					!	LOOP DIRECTIONS
					a_nm_mn					=	cmplx(0.0_dp,0.0_dp,dp)	
					do b = 1, 3
						a_nm_mn(:,b)		=	A_ka(:,n,m) 	* 	A_ka(b,m,n)
					end do
					a_nm_mn					=	dE_mn * a_nm_mn
					! --
					!
					!	LOOP HW
					d_broad(:)				= 	delta_broad(dE_mn, hw_lst(:),  i_eta_smr)
					a_dE_hw					=	cmplx(0.0_dp,0.0_dp,dp)
					do hw_idx = 1, n_hw
						a_dE_hw(:,:,hw_idx)	=	a_nm_mn(:,:) * d_broad(hw_idx)
					end do
					! --
					!
					!	LOOP FERMI LEVEL
					do ef_idx = 1, n_ef
						opt_herm(:,:,:,ef_idx)	=		opt_herm(:,:,:,ef_idx)		&
													- pi_dp * a_dE_hw(:,:,:) * df_mn(ef_idx)
					end do
					! --
				end if
			end do
		end do
		!$OMP END PARALLEL DO
		!
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
		!$OMP PARALLEL DO DEFAULT(	none )												&
		!$OMP PRIVATE(df_mn, dE_mn, b, a_nm_mn, a_dE_hw, hw_idx, ef_idx)			&
		!$OMP SHARED(n_wf,n_hw,n_ef, fd_distrib, en_k, A_ka, hw_lst, i_eta_smr )	&
		!$OMP REDUCTION(+: opt_aherm)													
		do m = 1, n_wf
			do n = 1, n_wf
				if(n/=m) then
					!
					!	dE	BANDS
					df_mn(:)	=	fd_distrib(:,m)		-	fd_distrib(:,n)
					dE_mn		=		en_k(m) 		- 		en_k(n)	
					! --
					!
					!	LOOP DIRECTIONS
					a_nm_mn	=	0.0_dp
					do b = 1, 3
						a_nm_mn(:,b)	=	A_ka(:,n,m) * A_ka(b,m,n)	
					end do
					a_nm_mn		=	a_nm_mn  * dE_mn
					! --
					!
					!	LOOP HW LASER
					a_dE_hw		=	cmplx(0.0_dp,0.0_dp,dp)
					do hw_idx = 1, n_hw
						a_dE_hw(:,:,hw_idx)	=	a_nm_mn(:,:)							&
												/	real(dE_mn - hw_lst(hw_idx) - i_eta_smr , dp)
					end do
					! --
					!
					!	LOOP FERMI LEVEL
					do ef_idx = 1, n_ef
						opt_aherm(:,:,:,ef_idx)	=	opt_aherm(:,:,:,ef_idx)	+ i_dp	* a_dE_hw(:,:,:)	* df_mn(ef_idx)
					end do
					! --
					!
				end if
			end do
		end do
		!$OMP END PARALLEL DO
		!
		return
	end function


















































end module kubo