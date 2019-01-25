module kubo
	!
	!	calculates some of the Kubo formulas available in postw90
	!	see wann guide chapter 12 for details
	!		http://www.wannier.org/doc/user_guide.pdf
	!
	!
	use constants,		only:	dp, i_dp, pi_dp
	use input_paras,	only:	kubo_tol
	use statistics,		only:	fd_stat
	implicit none


	private
	public			::			kubo_ahc_tens,		&
								kubo_ohc_tens,		&
								velo_ahc_tens,		&
								kubo_opt_tens	



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

	pure function velo_ahc_tens(en_k, v_kab, hw, eFermi, T_kelvin, i_eta_smr) result( o_ahc)
		!
		!	use Wanxiangs expression to calculate the ohc
		!
		!		VIA VELOCITIES
		real(dp),		intent(in)		::	en_k(:), hw, eFermi, T_kelvin
		complex(dp), 	intent(in)		::	v_kab(:,:,:), i_eta_smr
		complex(dp)						::	o_ahc(3,3), en_denom, delta_fd
		integer							::	n, l, j
		!
		o_ahc	=	0.0_dp
		!
		do n = 1, size(en_k)			
			do l = 1, size(en_k)
				!
				!
				if( l/=n ) then
					en_denom	= 	(	en_k(n) - en_k(l)			)**2 		-		 (		hw + i_eta_smr		)**2	
					!
					if(		abs(en_denom) > kubo_tol 		)  then
						delta_fd		=	fd_stat(en_k(n), eFermi, T_kelvin)	-	fd_stat(en_k(l), eFermi, T_kelvin)
						!
						do j = 1, 3
							o_ahc(:,j)	= 	o_ahc(:,j) 	+	 delta_fd * aimag(	v_kab(:,n,l) * v_kab(j,l,n)		) 	/		en_denom		
						end do
					end if
				end if
				!
				!
			end do
		end do
		!
		!
		return
	end function




	pure function kubo_ohc_tens(en_k, v_kab, hw, eFermi, T_kelvin, i_eta_smr)	result(z_ohc)
		!
		!		calculates wann guide (12.5) explicitly
		!	
		!	via VELOCITIES
		real(dp),		intent(in)		::	en_k(:), hw, eFermi, T_kelvin
		complex(dp),	intent(in)		::	v_kab(:,:,:), i_eta_smr
		complex(dp)						::	z_ohc(3,3)
		integer							::	n,	m, j
		real(dp)						::	dFdE_mn, dE_mn
		!
		z_ohc	=	cmplx(0.0_dp, 0.0_dp, dp)
		!
		do m = 1 , size(v_kab,3)
			do n = 1, size(v_kab,2)
				dE_mn 	=	en_k(m)	-	en_k(n)	
				dFdE_mn	=	fd_stat(en_k(m), 		eFermi, T_kelvin)	- fd_stat(en_k(n), 		eFermi, T_kelvin)
				!
				if( dE_mn > kubo_tol) then
					dFdE_mn	=	dFdE_mn	/	dE_mn
				end if
				!
				do j = 1, 3
					z_ohc(:,j)	=	z_ohc(:,j)	 + 	i_dp * 	cmplx(dFdE_mn,0.0_dp,dp)	* v_kab(:,n,m) * v_kab(j,m,n)	&
													/ 	(	cmplx(dE_mn - hw,0.0_dp,dp) - i_eta_smr	)
				end do
			end do
		end do
		!
		return
	end function




	subroutine kubo_opt_tens(hw, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka, opt_symm, opt_asymm)	
		!	see wann guide chapter 12
		!		eq. (12.14) & (12.15)
		real(dp),		intent(in)		::	hw, e_fermi, T_kelvin, en_k(:) 
		complex(dp),	intent(in)		::	i_eta_smr
		complex(dp), allocatable,	intent(in)		::	A_ka(:,:,:)
		complex(dp), 	intent(out)		::	opt_symm(3,3), 	opt_asymm(3,3)
		complex(dp)						::	opt_herm(3,3),	opt_aherm(3,3)
		!
		opt_symm	= 	0.0_dp
		opt_asymm	=	0.0_dp
		!
		if(			abs(hw)				< 1e-3_dp	&	!)	'[kubo_opt_tens]:	WARNING hw is very small: ',hw 
			.or.	abs(i_eta_smr)		< 1e-5_dp	&	!)	'[kubo_opt_tens]:	WARNING i_eta_smr is very small: ', i_eta_smr
			.or.	abs(T_kelvin)		< 1e-3_dp	&
		)then
			if(allocated(A_ka))	then
				write(*,*)	"[kubo_opt_tens]: 	WARNING opt_tens was set to zero; please make sure hw, i_eta_smr, T_kelvin are non zero"
			end if
		else
			if(allocated(A_ka)) then
				opt_herm	=	get_hermitian(hw, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)
				opt_aherm	=	get_anti_herm(hw, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)
				!
				opt_symm	=	real(	opt_herm	,	dp) 		+	i_dp * aimag(	opt_aherm	)		!	(12.14)		
				opt_asymm	=	real(	opt_aherm	,	dp)			+	i_dp * aimag(	opt_herm	)		!	(12.15)
			end if
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


	function get_hermitian(hw, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)	result(opt_herm)
		!	wann guide eq. (12.10)
		real(dp),		intent(in)		::	hw, e_fermi, T_kelvin, en_k(:)
		complex(dp),	intent(in)		::	i_eta_smr
		complex(dp),	intent(in)		::	A_ka(:,:,:)
		complex(dp)						::	opt_herm(3,3)
		integer							::	b, n, m 
		real(dp)						::	f_mn, dE_mn, pre_fact
		!
		opt_herm	=	0.0_dp
		!
		do m = 1, size(A_ka,3)
			do n = 1, size(A_ka,2)
				!
				f_mn	=	fd_stat(	en_k(m)		,e_fermi,T_kelvin)		-		fd_stat(	en_k(n)		,e_fermi,T_kelvin)
				dE_mn	=					en_k(m)								-						en_k(n) 
				!
				do b = 1, 3
					opt_herm(:,b)	=	opt_herm(:,b)	+	f_mn * dE_mn * 			A_ka(:,n,m) * A_ka(b,m,n)		* delta_broad( dE_mn - hw	, i_eta_smr)
				end do
				!
			end do
		end do
		!
		!
		pre_fact	=	-	pi_dp 
		opt_herm	=	cmplx(pre_fact,0.0_dp, dp)	* opt_herm
		!
		return
	end function



	function get_anti_herm(hw, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)	result(opt_aherm)
		!	wann guide eq. (12.11)
		real(dp),		intent(in)		::	hw, e_fermi, T_kelvin, en_k(:)
		complex(dp),	intent(in)		::	i_eta_smr	
		complex(dp),	intent(in)		::	A_ka(:,:,:)
		complex(dp)						::	opt_aherm(3,3)
		integer							::	b, n, m 
		real(dp)						::	f_mn, dE_mn 
		complex(dp)						::	pre_fact
		!
		opt_aherm	=	cmplx(0.0, 0.0, dp)
		!
		do m = 1, size(A_ka,3)
			do n = 1, size(A_ka,2)
				!
				f_mn	=		fd_stat(	en_k(m)		,e_fermi,T_kelvin)		-		fd_stat(	en_k(n)	 ,e_fermi, T_kelvin)
				dE_mn	=				en_k(m) - en_k(n)	
				!
				do b = 1, 3
					opt_aherm(:,b)	=	opt_aherm(:,b)	+	f_mn * 	real(		dE_mn / (dE_mn - hw - i_eta_smr)	,dp	) *			A_ka(:,n,m) * A_ka(b,m,n)		
				end do
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