module kubo
	!
	!	calculates some of the Kubo formulas available in postw90
	!	see wann guide chapter 12 for details
	!		http://www.wannier.org/doc/user_guide.pdf
	!
	!
	use constants,		only:	dp, i_dp, pi_dp, aUtoEv
	use input_paras,	only:	kubo_tol
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

	pure function kubo_ahc_tens(en_k, Om_kab, eFermi, T_kelvin) result( o_ahc)
		!	see wann guide chapter 12
		!		eq. (12.16)
		real(dp),		intent(in)		::	en_k(:), eFermi, T_kelvin
		complex(dp), allocatable,	intent(in)		::	Om_kab(:,:,:,:)
		real(dp)						::	o_ahc(3,3)
		integer							::	n
		!
		o_ahc	=	0.0_dp
		!
		if(allocated(Om_kab)) then
			do n = 1, size(en_k)
				o_ahc	=	o_ahc	-	real(Om_kab(:,:,n,n),dp)		*	fd_stat(en_k(n),	eFermi, 	T_kelvin) 
			end do
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

	pure function velo_ahc_tens(en_k, v_kab, hw_lst, eFermi, T_kelvin, i_eta_smr) result( o_ahc)
		!
		!	use Wanxiangs expression to calculate the ohc
		!
		!		VIA VELOCITIES
		real(dp),		intent(in)		::	en_k(:), hw_lst(:), eFermi, T_kelvin
		complex(dp), 	intent(in)		::	v_kab(:,:,:), i_eta_smr
		complex(dp)	,	allocatable		::	o_ahc(:,:,:)
		real(dp)						::	v_nm_mn(3,3), delta_fd, en_denom
		integer							::	n, l, j, hw
		!
		!
		allocate(	o_ahc(3,3,size(hw_lst))			)
		o_ahc	=	0.0_dp
		!
		do n = 1, size(en_k)
			do l = 1, size(en_k)
				if(l/=n) then
					!
					!
					en_denom		= 	(	en_k(n) - en_k(l)			)**2 	
					delta_fd		=	fd_stat(en_k(n), eFermi, T_kelvin)	-	fd_stat(en_k(l), eFermi, T_kelvin)
					!
					do j = 1, 3
						v_nm_mn(:,j)	=	aimag(	v_kab(:,n,l) * v_kab(j,l,n)		) 
					end do
					v_nm_mn			=	delta_fd * v_nm_mn
					!
					!
					do hw = 1, size(hw_lst,1)
						o_ahc(:,:,hw)	= 	o_ahc(:,:,hw) 	+	 v_nm_mn(:,:)	/		(	en_denom	- (		hw_lst(hw) + i_eta_smr		)**2	)	
					end do
					!
					!
				end if
			end do
		end do
		!
		!
		return
	end function




	pure function kubo_ohc_tens(en_k, v_kab, hw_lst, eFermi, T_kelvin, i_eta_smr)	result(z_ohc)
		!
		!		calculates wann guide (12.5) explicitly
		!	
		!	via VELOCITIES
		real(dp),		intent(in)		::	en_k(:), hw_lst(:), eFermi, T_kelvin
		complex(dp),	intent(in)		::	v_kab(:,:,:), i_eta_smr
		complex(dp), 	allocatable		::	z_ohc(:,:,:) 
		complex(dp)						::	v_nm_mn(3,3)
		integer							::	n,	m, j, hw
		real(dp)						::	dFdE_mn, dE_mn
		!
		allocate(	z_ohc(3,3,	size(hw_lst,1)		))
		z_ohc	=	cmplx(0.0_dp, 0.0_dp, dp)
		!
		do m = 1 , size(v_kab,3)
			do n = 1, size(v_kab,2)
				if(m/=n) then
					!
					!
					dE_mn 	=	en_k(m)	-	en_k(n)	
					dFdE_mn	=	fd_stat(en_k(m), 		eFermi, T_kelvin)	- fd_stat(en_k(n), 		eFermi, T_kelvin)
					dFdE_mn	=	dFdE_mn	/	dE_mn
					!
					do j = 1, 3
						v_nm_mn(:,j)	=	v_kab(:,n,m) * v_kab(j,m,n)
					end do
					v_nm_mn				=	v_nm_mn	*	dFdE_mn
					!
					!
					do hw = 1, size(hw_lst)
						z_ohc(:,:,hw)	=	z_ohc(:,:,hw)	 + 			i_dp * 	v_nm_mn(:,:)					&
																	/ 	(	dE_mn - hw_lst(hw) - i_eta_smr	)
					end do
					!
					!
				end if
			end do
		end do
		!
		return
	end function





	subroutine kubo_opt_tens(hw_lst, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka, opt_symm, opt_asymm)	
		!	see wann guide chapter 12
		!		eq. (12.14) & (12.15)
		!
		!	via CONNECTION
		real(dp),					intent(in)		::	hw_lst(:), e_fermi, T_kelvin, en_k(:) 
		complex(dp),				intent(in)		::	i_eta_smr
		complex(dp), allocatable,	intent(in)		::	A_ka(:,:,:)
		complex(dp), allocatable,	intent(out)		::	opt_symm(:,:,:), 	opt_asymm(:,:,:)
		complex(dp), allocatable					::	opt_herm(:,:,:),	opt_aherm(:,:,:)
		!
		allocate(	opt_symm(	3,3,size(hw_lst,1))		)
		allocate(	opt_asymm(	3,3,size(hw_lst,1))		)
		allocate(	opt_aherm(	3,3,size(hw_lst,1))		)
		allocate(	opt_herm(	3,3,size(hw_lst,1))		)
		opt_symm	= 	0.0_dp
		opt_asymm	=	0.0_dp
		!
		if(allocated(A_ka)) then
			opt_herm	=	get_hermitian(hw_lst, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)
			opt_aherm	=	get_anti_herm(hw_lst, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)
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



	function get_hermitian(hw_lst, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)	result(opt_herm)
		!	wann guide eq. (12.10)
		real(dp),		intent(in)		::	hw_lst(:), e_fermi, T_kelvin, en_k(:)
		complex(dp),	intent(in)		::	i_eta_smr
		complex(dp),	intent(in)		::	A_ka(:,:,:)
		complex(dp),	allocatable		::	opt_herm(:,:,:)
		complex(dp)						::	a_nm_mn(3,3)
		integer							::	n, m, b, hw 
		real(dp)						::	f_mn, dE_mn, pre_fact
		!
		allocate(opt_herm(3,3,size(hw_lst,1)))
		opt_herm	=	0.0_dp
		!
		do m = 1, size(A_ka,3)
			do n = 1, size(A_ka,2)
				if(n/=m) then
					!
					!
					f_mn	=	fd_stat(	en_k(m)		,e_fermi,T_kelvin)	-	fd_stat(	en_k(n)		,e_fermi,T_kelvin)
					dE_mn	=				en_k(m)							-				en_k(n) 
					!
					do b = 1, 3
						a_nm_mn(:,b)	=	A_ka(:,n,m) * A_ka(b,m,n)
					end do
					a_nm_mn				=	f_mn * dE_mn * a_nm_mn
					!
					!
					do hw = 1, size(hw_lst,1)
						opt_herm(:,:,hw)	=	opt_herm(:,:,hw)	+	a_nm_mn(:,:)		* delta_broad( dE_mn - hw_lst(hw)	, i_eta_smr)
					end do
					!
					!
				end if
			end do
		end do
		!
		!
		pre_fact	=	-	pi_dp 
		opt_herm	=	cmplx(pre_fact,0.0_dp, dp)	* opt_herm
		!
		return
	end function



	function get_anti_herm(hw_lst, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)	result(opt_aherm)
		!	wann guide eq. (12.11)
		real(dp),		intent(in)		::	hw_lst(:), e_fermi, T_kelvin, en_k(:)
		complex(dp),	intent(in)		::	i_eta_smr	
		complex(dp),	intent(in)		::	A_ka(:,:,:)
		complex(dp),	allocatable		::	opt_aherm(:,:,:)
		complex(dp)						::	a_nm_mn(3,3)
		integer							::	n, m, b, hw 
		real(dp)						::	f_mn, dE_mn 
		complex(dp)						::	pre_fact
		!
		allocate(	opt_aherm(3,3,size(hw_lst,1)))
		opt_aherm	=	cmplx(0.0, 0.0, dp)
		!
		do m = 1, size(A_ka,3)
			do n = 1, size(A_ka,2)
				!
				f_mn	=		fd_stat(	en_k(m)		,e_fermi,T_kelvin)		-		fd_stat(	en_k(n)	 ,e_fermi, T_kelvin)
				dE_mn	=				en_k(m) - en_k(n)	
				!
				do b = 1, 3
					a_nm_mn(:,b)	=	A_ka(:,n,m) * A_ka(b,m,n)	
				end do
				a_nm_mn		=	a_nm_mn	* f_mn * dE_mn
				!
				!
				do hw = 1, size(hw_lst,1)
					opt_aherm(:,:,hw)	=	opt_aherm(:,:,hw)	+				a_nm_mn(:,:)							&
																	/	real(dE_mn - hw_lst(hw) - i_eta_smr , dp)
				end do
				!
				!
			end do
		end do
		!
		!
		pre_fact	=	i_dp
		opt_aherm	=	pre_fact	*	opt_aherm
		!
		return
	end function


















































end module kubo