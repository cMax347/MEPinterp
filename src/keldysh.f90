module keldysh
	!

	!
	use omp_lib
	use input_paras,	only:	kubo_tol
	use constants,		only:	dp, i_dp, pi_dp, aUtoEv
	use statistics,		only:	fd_stat


	implicit none


	private
	public			::		keldysh_scnd_photoC

	real(dp),	parameter		::	pi_half	= pi_dp / 2.0_dp, &
									band_tol	=	1e-5_dp
contains


!
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				***		PUBLIC KELDYSH FUNCTIONS 		***																							|
!				->	2nd order photocurrents																											|
!																																					|
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
function keldysh_scnd_photoC(	en_k, V_ka, hw_lst, smr_lst, ef_lst)			result(phi)
	!		see	 Freimuth et al., PRB 94, 144432 (2016)  
	!		 EQ.(2)
	real(dp),		intent(in)		::	en_k(:), hw_lst(:), smr_lst(:), ef_lst(:)
	complex(dp),	intent(in)		::	V_ka(:,:,:)
	complex(dp),	allocatable		::	phi(:,:,:,:,:,:)
	complex(dp)						::	v_ijk(3,3,3), v_ikj(3,3,3)
	integer							::	a,b,c, j,k, hw, smr, ef
	!
	allocate(	phi(	3,3,3,	size(hw_lst), size(smr_lst), size(ef_lst)		))
	phi	= cmplx(0.0_dp,0.0_dp,dp)
	!
	!	TODO:	OMP REDUCTION OVER "phi"
	do c = 1, size(en_k,1)
		do b = 1, size(en_k,1)
			do a = 1, size(en_k,1)
				!
				!	GET VELO NOMINATOR
				do k = 1, 3
					do j = 1,3
						v_ijk(:,j,k)	=		V_ka(:,c,a) * V_ka(j,a,b) * V_ka(k,b,c)	
						!
						v_ikj(:,k,j)	=		V_ka(:,c,a) * V_ka(k,a,b) * V_ka(j,b,c)	
					end do
				end do
				!
				!	GET ENERGY DENOM
				do ef = 1, size(ef_lst)
					do smr = 1, size(smr_lst)
						do hw = 1, size(hw_lst) 
							!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
							phi(:,:,:,hw,smr,ef)=phi(:,:,:,hw,smr,ef)																				&
								+(			v_ijk(:,:,:)																							&
										* en_int_RRR(	en_k(a),	en_k(b)-hw_lst(hw)		,en_k(c),		ef_lst(ef)			, smr_lst(smr))		&
									!~~~~~~
									-		v_ijk(:,:,:)																							&
										* en_int_RRA(	en_k(a),	en_k(b)-hw_lst(hw)		,en_k(c),		ef_lst(ef)			, smr_lst(smr))		&
									!~~~~~~
									+		v_ikj(:,:,:)																							&
										* en_int_RRR(	en_k(a),	en_k(b)+hw_lst(hw)		,en_k(c),		ef_lst(ef)			, smr_lst(smr))		&
									!~~~~~~
									-		v_ikj(:,:,:)																							&
										* en_int_RRA(	en_k(a),	en_k(b)+hw_lst(hw)		,en_k(c),		ef_lst(ef)			, smr_lst(smr))		&
									!~~~~~~
									+		v_ijk(:,:,:)																							&
										* en_int_RRA(	en_k(a),	en_k(b)-hw_lst(hw)		,en_k(c),	ef_lst(ef)-hw_lst(hw)	, smr_lst(smr))		&
									!~~~~~~
									+		v_ikj(:,:,:)																							&
										* en_int_RRA(	en_k(a),	en_k(b)+hw_lst(hw)		,en_k(c),	ef_lst(ef)+hw_lst(hw)	, smr_lst(smr))		&
								) / hw_lst(hw)**2	
							!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
						end do
					end do
				end do
				!
				!
			end do
		end do
	end do
	!
	return
end function




complex(dp) function keldysh_torques()
	!
	!		ToDo
	!
	keldysh_torques	=	0.0_dp
	!
	return
end function
!~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






!
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				***		PRIVATE HELPERS - ANALYTIC ENERGY INTEGRATION		***																		|
!																																					|
!																																					|
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				***		"G_R G_R G_R" 	INTEGRALS	***										|
!																							|
!																							|
complex(dp) pure function en_int_RRR(E1, E2, E3, E4,  smr)
	!		see	 Freimuth et al., PRB 94, 144432 (2016)  
	!		 EQ.(B7)	/	EQ.(B9)		/ EQ.(B11)
	!		i.e.	I1(E1,E2,E3,E4)
	real(dp),		intent(in)		::	E1, E2, E3, E4, smr
	real(dp)						::	smrsmr
	!
	!
	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	!			BAND DEGENERACY
	!				use EQ.(B9)/(B11)
	!	with identity
	!	->	EQ.(B9):		I1(E1,E1,E3,E4)	=	I1(E1,E3,E1,E4)	=	I1(E3,E1,E1,E4)
	!	->	EQ.(B11):			I1(E1,E1,E1,E4)	
	!
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	if(					(	abs(E1-E2)	<	band_tol ) 			&	!	E1==E2
				.and.	(	abs(E2-E3)	>	band_tol )			&	!	E2/=E3
	)then
		!	Eq.(B9):
			en_int_RRR	=	en_int_DGN_RRR(	E1, E3, E4, smr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E1-E3)	<	band_tol ) 			&	!	E1==E3
				.and.	(	abs(E2-E3)	>	band_tol )			&	!	E2/=E3)then
	)then
		!	Eq.(B9):
			en_int_RRR	=	en_int_DGN_RRR( E1, E2, E4, smr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E2-E3)	<	band_tol ) 			&	!	E2==E3
				.and.	(	abs(E1-E3)	>	band_tol )			&	!	E1/=E3)then)then
	)then
		!	Eq.(B9):
			en_int_RRR	=	en_int_DGN_RRR(E2, E1, E4, smr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E2-E3)	<	band_tol )			&	!	E1==E2==E3																											
				.and.	(	abs(E1-E3)	<	band_tol )			&	!	E1==E2==E3																										
	)then
		!	Eq.(B11):																												!
			en_int_RRR	=	-1.0_dp	/	(2.0_dp*( E4-E1	+ i_dp*smr	)**2)																			!
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	else
		smrsmr	= smr**2
		!	EQ.(B7):	 I1(E1,E2,E3,E4)
		!		
		en_int_RRR	=	&
						log(	1.0_dp	+	(E1 - E4)**2	/smrsmr	)	/	(	2.0_dp	* (E1-E2) * (E1-E3)	)				&
					+	log(	1.0_dp	+	(E2 - E4)**2	/smrsmr	)	/	(	2.0_dp	* (E2-E3) * (E2-E1)	)				&
					+	log(	1.0_dp	+	(E3 - E4)**2	/smrsmr	)	/	(	2.0_dp	* (E3-E1) * (E3-E2)	)				&
					!
					+	atan(	(	E4	-	E1	)			/  smr	)	/	(	i_dp 	* (E1-E2) * (E1-E3) )				&
					+	atan(	(	E4	-	E2	)			/  smr	)	/	(	i_dp 	* (E2-E3) * (E2-E1) )				&
					+	atan(	(	E4	-	E3	)			/  smr	)	/	(	i_dp 	* (E3-E1) * (E3-E2) )											
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	end if
	!
	return
end function 
!~
!~
complex(dp) pure function en_int_DGN_RRR(E1, E3, E4, smr)
	real(dp),		intent(in)		::	E1,E3,E4,smr
	!		see	 Freimuth et al., PRB 94, 144432 (2016)  
	!		 EQ.(B9)
	!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	en_int_DGN_RRR	=	&
						log( 	( smr**2 + (E3 - E4)**2	)	/	( smr**2 + (E1 - E4)**2 ))		/	(	2.0_dp 	* (E1-E3)**2	)	&	
						+	(		atan((E4-E3)/smr)		-		atan((E4-E1)/smr)	)		/	(	i_dp 	* (E1-E3)**2 	)	&
						+	1.0_dp	/	(	(E3-E1)*(E4-E1+i_dp*smr)	)

	return
end function
!~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~





!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				***		"G_R G_R G_A" 	INTEGRALS	***										|
!																							|
!																							|
complex(dp) pure function en_int_RRA(E1, E2, E3, E4,  smr)
	!		see	 Freimuth et al., PRB 94, 144432 (2016) 
	!		 EQ.(B8)	/	EQ.(B10)	
	!		i.e.	I2(E1,E2,E3,E4)
	real(dp),		intent(in)		::	E1, E2, E3, E4, smr
	complex(dp)						::	i_2smr,smrsmr
	!	
	i_2smr		=	2.0_dp * i_dp * smr	
	smrsmr		=	smr**2	
	en_int_RRA	=	cmplx(0.0_dp,0.0_dp,dp)						
	!
	if(					(	abs(E1-E2)	<	band_tol )&
				.and.	(	abs(E2-E3)	>	band_tol )&		
	)then
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
		!	Eq.(B10):			I2(E1,E1,E3,E4)																										!
			en_int_RRA	=	en_int_RRA																						&
						 + log( 	(smrsmr+(E3-E4)**2)/	(smrsmr+(E1-E4)**2))/		(	2.0_dp* (E1-E3-i_2smr)**2	)	&						!
						 + i_dp * (		pi_half 	+	atan((E4-E3)/smr)	)	/		(			(E1-E3-i_2smr)**2	)	&						!
						 + i_dp * (		pi_half		+	atan((E4-E1)/smr)	)	/		(			(E1-E3-i_2smr)**2	)	&						!
						 +  1.0_dp												/		((E3-E1+i_2smr)*(E4-E1+i_dp*smr))							!
		!___________________________________________________________________________________________________________________________________________!
	else
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
		!	Eq.(B8):				I2(E1,E2,E3,E4)																									!
																																					!																																	!
		en_int_RRA	=	en_int_RRA																												 &	!
					+	log(	1.0_dp	+	(E1 - E4)**2 	/smrsmr		)		/	(	2.0_dp *	 (E1-E2) 		*	( E1-E3 -i_2smr)		)&	!		
					+	log(	1.0_dp	+	(E2 - E4)**2	/smrsmr		)		/	(	2.0_dp * ( E2-E3 -i_2smr) 	* 		(E2-E1)				)&	!
					+	log(	1.0_dp	+	(E3 - E4)**2	/smrsmr		)		/	(	2.0_dp * ( E3-E1 +i_2smr) 	*	( E3-E2 +i_2smr)		)&	!
					!-----																															!
					+	i_dp * (	pi_half		+	atan((E4-E1)/smr)	)		/	(			(E1-E2) 			*	( E3-E1 + i_2smr)		)&	!
					+	i_dp * (	pi_half		+	atan((E4-E2)/smr)	)		/	(		( E3-E2 + i_2smr)		*		(E2-E1)				)&	!
					+	i_dp * (	pi_half		+	atan((E4-E3)/smr)	)		/	(		( E3-E1 + i_2smr)		*	( E3-E2 + i_2smr)		)	!					
		!___________________________________________________________________________________________________________________________________________!
	end if
	!
	return
end function 
!~
!~













end module keldysh










