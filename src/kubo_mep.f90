module kubo_mep
	!
	use constants,		only:	dp, i_dp, pi_dp
	use matrix_math,	only:	my_Levi_Civita, convert_tens_to_vect
	use input_paras,	only:	kubo_tol
	use statistics,		only:	fd_stat
	implicit none


	private
	public			::			kubo_mep_CS,		&
								kubo_mep_LC,		&
								kubo_mep_IC
							



contains

	pure function kubo_mep_CS(eFermi, T_kelvin, en_k,  A_ka, Om_kab) result(cs_tens)
		real(dp),						intent(in)		::	eFermi, T_kelvin, en_k(:)
		complex(dp),	allocatable, 	intent(in)		::	A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp)										::	cs_tens(3,3)
		complex(dp)										::	cs_scal, om_vect(3)
		integer											::	n0, i
		!
		cs_tens	=	0.0_dp
		if(	allocated(A_ka)		.and.		allocated(Om_kab)		)then
			cs_scal	= 0.0_dp
			!
			!	sum over valence
			do n0 = 1, size(Om_kab,3)
				call convert_tens_to_vect(	Om_kab(:,:,n0,n0),		om_vect(:)	)
				cs_scal	= cs_scal +	 0.5_dp  	* 		fd_stat(en_k(n0),eFermi, T_kelvin)  	*	dreal(	 dot_product(	A_ka(:,n0,n0)	,	om_vect(:)	)			)
			end do
			!
			do i = 1, 3
				cs_tens(i,i)	= real(cs_scal,dp)
			end do
		end if
		!
		return
	end function


	
	function kubo_mep_LC(eFermi, T_kelvin, velo, En, neglected) result(F3)
		!	LOCAL KUBO CONTRIBUTION
		!
		real(dp),			intent(in)		::	eFermi, T_kelvin, En(:)
		complex(dp),		intent(in)		::	velo(:,:,:)
		integer,			intent(out)		::	neglected
		real(dp)							::	F3(3,3)
		integer								::	n0, n, tot, 		&
												i, j, k, l	
		real(dp) 							::	pre_fact, velo_nom, en_term, en_denom, fermi_dirac
		!
		neglected	= 0
		tot 		= 0
		!
		F3	= 0.0_dp
		do n0 = 1, size(velo,2)
			!
			fermi_dirac	= fd_stat(en(n0), eFermi, T_kelvin)
			!
			!MIXING
			do n = 1, size(velo,2)
			 	if( n/= n0 )	then
			 		en_denom	=	(	en(n0) - en(n)	)**3
			 		en_term		=	fermi_dirac	/	en_denom		
			 		tot			= 	tot + 1
			 		!
					!TRIPLE PRODUCT			 		
			 		do j = 1, 3
			 			do i = 1, 3
			 				do l = 1, 3
			 					do k = 1, 3
			 						pre_fact	= 	real(my_Levi_Civita(j,k,l),dp)  
			 						!
			 						velo_nom	=	real(		velo(i,n0,n) * velo(k,n,n0) * velo(l,n0,n0)		, dp)
			 						!		
			 						!
			 						F3(i,j)		= 	F3(i,j) 	+  		pre_fact * 	velo_nom  * en_term		
			 					end do
			 				end do
			 			end do
			 		end do
			 	end if
			 	!
			 	!
			end do
		end do
		!
		return
	end function


	function kubo_mep_IC(eFermi, T_kelvin, velo, En, neglected) result(F2)
		!	ITINERANT KUBO CONTRIBUTION
		!
		complex(dp),		intent(in)		::	velo(:,:,:)
		real(dp),			intent(in)		::	eFermi, T_kelvin, En(:)
		integer,			intent(out)		::	neglected
		real(dp)							::	F2(3,3)
		integer								::	n0, m, n, 		&
												i, j, k, l,		&
												tot
		real(dp) 							::	pre_fact, velo_nom, en_term, en_denom, fermi_dirac
		!
		neglected	= 0
		tot 		= 0
		!
		F2	= 0.0_dp
		do n0 = 1, size(velo,2)
			!
			fermi_dirac	= fd_stat(en(n0), eFermi, T_kelvin)
			!
			!MIXING
			do m = 1, size(velo,2)
				do n = 1, size(velo,2)
				 	if( n/= n0 .and. m/=n0 )	then
				 		en_denom	=	(	en(n0) - en(n)	)**2		 * 		(	en(n0) - en(m)	)  
				 		en_term		=	fermi_dirac		/ 	en_denom
				 		tot 		= 	tot + 1
				 		!
				 		!TRIPLE PRODUCT
				 		do j = 1, 3
				 			do i = 1, 3
				 				do l = 1, 3
				 					do k = 1, 3
				 						pre_fact	=	- real(my_Levi_Civita(j,k,l),dp) 
				 						!
				 						velo_nom	=	real(	velo(i,n0,n) * velo(k,n,m) * velo(l,m,n0)		, dp)
				 						!
				 						!
				 						F2(i,j)		= 	F2(i,j)		+	pre_fact  * velo_nom	* en_term
				 					end do
				 				end do
				 			end do
				 		end do
				 	end if
				 	!
				 	!
				end do
			end do
		end do
		!
		return
	end function




end module kubo_mep