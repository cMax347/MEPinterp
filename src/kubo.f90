module kubo
	!
	!	calculates some of the Kubo formulas available in postw90
	!	see wann guide chapter 12 for details
	!		http://www.wannier.org/doc/user_guide.pdf
	!
	!
	use constants,		only:	dp, i_dp, pi_dp
	use statistics,		only:	fd_stat
	implicit none


	private
	public			::			kubo_ahc_tens,		&
								velo_ahc_tens,		&
								kubo_opt_tens	



contains





!public
	pure function kubo_ahc_tens(en_k, Om_kab, eFermi, T_kelvin, unit_vol) result( o_ahc)
		!	see wann guide chapter 12
		!		eq. (12.16)
		real(dp),		intent(in)		::	en_k(:), eFermi, T_kelvin, unit_vol
		complex(dp),	intent(in)		::	Om_kab(:,:,:,:)
		real(dp)						::	o_ahc(3,3)
		integer							::	n
		!
		o_ahc	=	0.0_dp
		!
		do n = 1, size(en_k)
			o_ahc	=	o_ahc	-	real(Om_kab(:,:,n,n),dp)		*	fd_stat(en_k(n),	eFermi, 	T_kelvin) /	unit_vol
		end do
		!
		!
		return
	end function

	pure function velo_ahc_tens(en_k, v_kab, eFermi, T_kelvin, unit_vol) result( o_ahc)
		real(dp),		intent(in)		::	en_k(:), eFermi, T_kelvin, unit_vol
		complex(dp), 	intent(in)		::	v_kab(:,:,:)
		real(dp)						::	o_ahc(3,3), en_denom	
		real(dp)						::	Om_ab(3,3)
		integer							::	n, l, j
		!
		o_ahc	=	0.0_dp
		!
		do n = 1, size(en_k)
			Om_ab	= 	0.0_dp
			!
			!	get curvature of band n
			do l = 1, size(en_k)
				en_denom	= en_k(n)-en_k(l)	
				!
				!
				if(		(l /= n) .and.	abs(en_denom) > 1e-6_dp 		)  then
					do j = 1, 3
						Om_ab(:,j)	= Om_ab(:,j) 	-2.0_dp * aimag(	v_kab(:,n,l) * v_kab(j,l,n)		) 	/	( 	en_denom	)**2 		
					end do
				end if
				!
			end do
			!
			!
			o_ahc	=	o_ahc	+ 	Om_ab 	*		fd_stat(en_k(n),	eFermi, T_kelvin)	/	unit_vol
		end do
		!
		!
		return
	end function




	subroutine kubo_opt_tens(hw, vol, e_fermi, T_kelvin, i_eta_smr , en_k, A_ka, opt_symm, opt_asymm)	
		!	see wann guide chapter 12
		!		eq. (12.14) & (12.15)
		real(dp),		intent(in)		::	hw, vol, e_fermi, T_kelvin, en_k(:) 
		complex(dp),	intent(in)		::	i_eta_smr
		complex(dp),	intent(in)		::	A_ka(:,:,:)
		complex(dp), 	intent(out)		::	opt_symm(3,3), 	opt_asymm(3,3)
		complex(dp)						::	opt_herm(3,3),	opt_aherm(3,3)
		!
		opt_herm	=	get_hermitian(hw, vol, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)
		opt_aherm	=	get_anti_herm(hw, vol, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)
		!
		!
		opt_symm	=	real(	opt_herm	,	dp) 		+	i_dp * aimag(	opt_aherm	)		!	(12.14)		
		opt_asymm	=	real(	opt_aherm	,	dp)			+	i_dp * aimag(	opt_herm	)		!	(12.15)
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
		return
	end function


	function get_hermitian(hw, unit_vol, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)	result(opt_herm)
		!	wann guide eq. (12.10)
		real(dp),		intent(in)		::	hw, unit_vol, e_fermi, T_kelvin, en_k(:)
		complex(dp),	intent(in)		::	i_eta_smr
		complex(dp),	intent(in)		::	A_ka(:,:,:)
		complex(dp)						::	opt_herm(3,3)
		integer							::	a, b, n, m 
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
					do a = 1, 3
						opt_herm(a,b)	=	opt_herm(a,b)	+	f_mn * dE_mn * 			A_ka(a,n,m) * A_ka(b,m,n)		* delta_broad( dE_mn - hw	, i_eta_smr)
					end do
				end do
				!
			end do
		end do
		!
		!
		pre_fact	=	-	pi_dp / unit_vol
		opt_herm	=	cmplx(pre_fact,dp)	* opt_herm
		!
		return
	end function



	function get_anti_herm(hw, unit_vol, e_fermi, T_kelvin, i_eta_smr, en_k, A_ka)	result(opt_aherm)
		!	wann guide eq. (12.11)
		real(dp),		intent(in)		::	hw, unit_vol, e_fermi, T_kelvin, en_k(:)
		complex(dp),	intent(in)		::	i_eta_smr	
		complex(dp),	intent(in)		::	A_ka(:,:,:)
		complex(dp)						::	opt_aherm(3,3)
		integer							::	a, b, n, m 
		real(dp)						::	f_mn, dE_mn 
		complex(dp)						::	pre_fact
		!
		opt_aherm	=	cmplx(0.0, 0.0, dp)
		!
		do m = 1, size(A_ka,3)
			do n = 1, size(A_ka,2)
				!
				f_mn	=		fd_stat(	en_k(m)		,e_fermi,T_kelvin)		-		fd_stat(	en_k(n)		,e_fermi,T_kelvin)
				dE_mn	=		real(		(		en_k(m)-en_k(n)		)		/	(	en_k(m)-en_k(n)-hw-i_eta_smr	)		,dp)
				!
				do b = 1, 3
					do a = 1,3 
						opt_aherm(a,b)	=	opt_aherm(a,b)	+	f_mn * 	dE_mn *			A_ka(a,n,m) * A_ka(b,m,n)		
					end do
				end do
				!
			end do
		end do
		!
		!
		pre_fact	=	i_dp /	unit_vol
		opt_aherm	=	pre_fact	*	opt_aherm
		!
		return
	end function


















































end module kubo