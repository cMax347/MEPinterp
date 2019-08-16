module keldysh
	!

	!
	use omp_lib
	use input_paras,	only:	kubo_tol
	use constants,		only:	dp, i_dp, pi_dp
	use statistics,		only:	fd_stat


	implicit none


	private
	public			::		keldysh_scnd_photoC

	real(dp),	parameter		::	pi_half	= pi_dp / 2.0_dp
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
									+		v_ikj(:,:,:)																					&
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
							!write(*,'(a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a,e12.4,a)')	"E1,E2-hw,E2+hw,E3,E4,E4-hw,E4+hw:",	&
							!													en_k(a)," ",en_k(b)-hw_lst(hw)," ",en_k(b)+hw_lst(hw),&
							!												" ",en_k(c)," ",ef_lst(ef)," ", ef_lst(ef)-hw_lst(hw)," ",ef_lst(ef)+hw_lst(hw)
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
complex(dp) function en_int_RRR(E1, E2, E3, E4,  smr)
	!		see	 Freimuth et al., PRB 94, 144432 (2016)  
	!		 EQ.(B7)	/	EQ.(B9)		/ EQ.(B11)
	!		i.e.	I1(E1,E2,E3,E4)
	real(dp),		intent(in)		::	E1, E2, E3, E4, smr
	real(dp)						::	smrsmr
	real(dp)						::	denom_1213, denom_2321, denom_3132, dE_14, dE_24, dE_34
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
	if(					(	abs(E1-E2)	<	kubo_tol ) 			&	!	E1==E2
				.and.	(	abs(E1-E3)	>	kubo_tol )			&	!	E1/=E3
	)then
		!	Eq.(B9):
			en_int_RRR	=	en_int_DGN_RRR(	E1, E3, E4, smr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E1-E3)	<	kubo_tol ) 			&	!	E1==E3
				.and.	(	abs(E2-E1)	>	kubo_tol )			&	!	E1/=E2)then
	)then
		!	Eq.(B9):
			en_int_RRR	=	en_int_DGN_RRR( E2, E1, E4, smr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E2-E3)	<	kubo_tol ) 			&	!	E2==E3
				.and.	(	abs(E1-E2)	>	kubo_tol )			&	!	E1/=E2)then)then
	)then
		!	Eq.(B9):
			en_int_RRR	=	en_int_DGN_RRR(E2, E1, E4, smr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E2-E3)	<	kubo_tol )			&	!	E1==E2==E3																											
				.and.	(	abs(E1-E3)	<	kubo_tol )			&	!	E1==E2==E3																										
	)then
		!	Eq.(B11):																												!
			en_int_RRR	=	-1.0_dp	/	(2.0_dp*( E4-E1	+ i_dp*smr	)**2)	
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	else
		smrsmr	= smr**2
		!	EQ.(B7):	 I1(E1,E2,E3,E4)
		
		dE_14		=	E4 -E1 
		dE_24		=	E4 -E2 
		dE_34		=	E4 -E3 
		!
		denom_1213	=	(E1-E2) * (E1-E3)
		denom_2321	=	(E2-E3) * (E2-E1)
		denom_3132	=	(E3-E1) * (E3-E2)
		!
		!
		en_int_RRR	=	&
						log(	1.0_dp	+	dE_14**2	/smrsmr	)	/	(	2.0_dp	* denom_1213	)		&
					+	log(	1.0_dp	+	dE_24**2	/smrsmr	)	/	(	2.0_dp	* denom_2321	)		&
					+	log(	1.0_dp	+	dE_34**2	/smrsmr	)	/	(	2.0_dp	* denom_3132	)		&
					!
					+	atan(	dE_14			/  smr	)	/	(	i_dp 	* denom_1213	)				&
					+	atan(	dE_24			/  smr	)	/	(	i_dp 	* denom_2321	)				&
					+	atan(	dE_34			/  smr	)	/	(	i_dp 	* denom_3132	)											
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	end if
	!
	!
	return
end function 
!~
!~
complex(dp) pure  function en_int_DGN_RRR(E1, E3, E4, smr)
	real(dp),		intent(in)		::	E1,E3,E4,smr
	real(dp)						::	dE_13, denom_14, re_RRR, im_RRR
	!		see	 Freimuth et al., PRB 94, 144432 (2016)  
	!		 EQ.(B9)
	!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	dE_13		=	E3-E1
	denom_14	=	(E4-E1)**2 + smr**2
	!
	re_RRR		=	0.5_dp * 	log( 	( smr**2 + (E3 - E4)**2	)	/	( smr**2 + (E1 - E4)**2 ))			&	
					- 	dE_13 *(E4-E1)	/ denom_14							
	!
	im_RRR		=	-atan((E4-E3)/smr)	&
					+atan((E4-E1)/smr)	& 
					-smr * dE_13 / denom_14	
	!
	re_RRR		=	re_RRR	/	dE_13**2
	im_RRR		=	im_RRR	/	dE_13**2
	!
	en_int_DGN_RRR	=	cmplx(	re_RRR,	im_RRR,		dp)
	!
	!^^^^^^^^^^^^^^^^^^^^^^^^
	!en_int_DGN_RRR	=	(		0.5_dp * 	log( 	( smr**2 + (E3 - E4)**2	)	/	( smr**2 + (E1 - E4)**2 ))			&	
	!						- 	dE_13 *(E4-E1)	/ denom_14																&
	!					)/ dE_13**2																						&
	!					-	i_dp   *(	atan((E4-E3)/smr)	-	atan((E4-E1)/smr) + smr* dE_13 / denom_14 )	/ dE_13**2		
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
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
	!
	if(					(	abs(E1-E2)	<	kubo_tol )&
				.and.	(	abs(E2-E3)	>	kubo_tol )&		
	)then
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
		!	Eq.(B10):			I2(E1,E1,E3,E4)																										!
			en_int_RRA	=	&																														!
						 + log( 	(smrsmr+(E3-E4)**2)/	(smrsmr+(E1-E4)**2))/		(	2.0_dp* (E1-E3-i_2smr)**2	)	&						!
						 + i_dp * (		pi_half 	+	atan((E4-E3)/smr)	)	/		(			(E1-E3-i_2smr)**2	)	&						!
						 + i_dp * (		pi_half		+	atan((E4-E1)/smr)	)	/		(			(E1-E3-i_2smr)**2	)	&						!
						 +  1.0_dp												/		((E3-E1+i_2smr)*(E4-E1+i_dp*smr))							!
		!___________________________________________________________________________________________________________________________________________!
		!en_int_RRA	=	cmplx(0.0_dp,0.0_dp,dp)		
	else
		!
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
		!	Eq.(B8):				I2(E1,E2,E3,E4)																									!
																																					!																																	!
		en_int_RRA	=	&																															!
					+	log(	1.0_dp	+	(E1 - E4)**2 	/smrsmr		)		/	(	2.0_dp *	 (E1-E2) 		*	( E1-E3 -i_2smr)		)&	!		
					+	log(	1.0_dp	+	(E2 - E4)**2	/smrsmr		)		/	(	2.0_dp * ( E2-E3 -i_2smr) 	* 		(E2-E1)				)&	!
					+	log(	1.0_dp	+	(E3 - E4)**2	/smrsmr		)		/	(	2.0_dp * ( E3-E1 +i_2smr) 	*	( E3-E2 +i_2smr)		)&	!
					!-----																															!
					+	i_dp * (	pi_half		+	atan((E4-E1)/smr)	)		/	(			(E1-E2) 			*	( E3-E1 + i_2smr)		)&	!
					+	i_dp * (	pi_half		+	atan((E4-E2)/smr)	)		/	(		( E3-E2 + i_2smr)		*		(E2-E1)				)&	!
					+	i_dp * (	pi_half		+	atan((E4-E3)/smr)	)		/	(		( E3-E1 + i_2smr)		*	( E3-E2 + i_2smr)		)	!					
		!___________________________________________________________________________________________________________________________________________!
		!en_int_RRA	=	cmplx(0.0_dp,0.0_dp,dp)		
	end if
	!
	!en_int_RRA	=	cmplx(0.0_dp,0.0_dp,dp)				
	!
	return
end function 
!~
!~













end module keldysh










