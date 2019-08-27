module green_en_int
	use omp_lib
	use input_paras,	only:	kubo_tol, E_int_MIN, N_energy_int
	use constants,		only:	dp, i_dp, pi_dp, aUtoEv	
	use statistics,		only:	fd_stat
	use matrix_math,	only:	mat_tr
	!use m_integrals_rrr,only:	integrals_rrr
	!use m_integrals_rra,only:	integrals_rra
	!
	implicit none
	!
	private
	public 		::				en_int_RRR,				&
								en_int_RRA
	!
	real(dp),	parameter	::	pi_half	=	pi_dp * 0.5_dp, &
								err_tol	=	1e-1_dp
	logical,	parameter	::	verbose_RRR	=	.True.,		&
								verbose_RRA	=	.False.
	!
contains




	complex(dp) pure function en_int_RRR(E1, E2, E3, E4, smr)
		real(dp),		intent(in)		::	E1, E2, E3, E4, smr
		complex							:: frank_RRR, my_RRR
		!
		!	call my routine
		my_RRR	=	my_en_int_RRR(E1,E2,E3, E4, smr)	
		!
		!!	call franks routine
		!call integrals_rrr(real(E1),real(E2),real(E3), real(E4),real(smr),	frank_RRR)
		!!	compare and return franks
		!if(abs(my_RRR-frank_RRR)> err_tol	.and. verbose_RRR ) then
		!	!write(*,'(a,e10.3,a,e10.3,a,f6.2,a,f6.2,a,f6.2,a,f6.2)')					&	
		!	write(*,'(a,e10.3,a,e12.6)')&
		!			"[en_int_RRR]: 	WARNING    abs error=",abs(my_RRR-frank_RRR),"	rel error=", abs(my_RRR-frank_RRR)/frank_RRR
		!			!" at dE_13="!,E3-E1," dE_14=",E4-E1," dE_12=",E2-E1," smr=",smr
		!end if 
		en_int_RRR 	=	my_RRR!frank_RRR
		!
		return
	end function


	complex(dp) pure function en_int_RRA(E1, E2, E3, E4, smr)
		real(dp),		intent(in)		::	E1, E2, E3, E4, smr
		complex							:: frank_RRA, my_RRA
		!
		!	call my routine
		my_RRA	=	my_en_int_RRA(E1,E2,E3, E4, smr)
		!!	call franks routine
		!call integrals_rra(real(E1),real(E2),real(E3), real(E4), real(smr),	frank_RRA)
		!!	compare and return franks
		!if(abs(my_RRA-frank_RRA)> err_tol .and. verbose_RRA )&	
		!	write(*,'(a,e10.3,a,f6.2,a,f6.2,a,f6.2,a,f6.2)')&	
		!			"[en_int_RRA]: 	WARNING   error=",abs(my_RRA-frank_RRA),"		 @ dE_13=",E3-E1," dE_14=",E4-E1," dE_12=",E2-E1," smr=",smr
		en_int_RRA 	=	my_RRA!frank_RRA
		!
		return
	end function



!
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				***		MY IMPLEMENTATION OF ANALYTIC INTEGRALS		***																				|
!																																					|
!																																					|
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				***		"G_R G_R G_R" 	INTEGRALS	***										|
!																							|
!																							|



	!
	complex(dp) pure function my_en_int_RRR(E1, E2, E3, E4,  smr)
	!		see	 Freimuth et al., PRB 94, 144432 (2016)  
	!		 EQ.(B7)	/	EQ.(B9)		/ EQ.(B11)
	!		i.e.	I1(E1,E2,E3,E4)
	real(dp),		intent(in)		::	E1, E2, E3, E4, smr
	real(dp)						::	smrsmr
	real(dp)						::	denom_1213, denom_2321, denom_3132, dE_14, dE_24, dE_34
	!
	!
	dE_14		=	E4 -E1
	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	!			BAND DEGENERACY
	!				use EQ.(B9)/(B11)
	!	with identity
	!	->	EQ.(B9):		I1(E1,E1,E3,E4)	=	I1(E1,E3,E1,E4)	=	I1(E3,E1,E1,E4)
	!	->	EQ.(B11):			I1(E1,E1,E1,E4)	
	!
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	if (				(	abs(E1-E2)	<	kubo_tol )			&
				.and.	(	abs(E1-E3)	<	kubo_tol )			&
				.and.	(	abs(E2-E3)	<	kubo_tol )			&
	)then
		!	Eq.(B11):																												!
		my_en_int_RRR	=	-1.0_dp	/	(2.0_dp*( dE_14	+ i_dp*smr	)**2)	
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E1-E2)	<	kubo_tol ) 			&	!	E1==E2
				.and.	(	abs(E1-E3)	>	kubo_tol )			&	!	E1/=E3
	)then
		!	Eq.(B9):
			my_en_int_RRR	=	en_int_DGN_RRR(	E1, E3, E4, smr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E1-E3)	<	kubo_tol ) 			&	!	E1==E3
				.and.	(	abs(E2-E1)	>	kubo_tol )			&	!	E1/=E2)then
	)then
		!	Eq.(B9):
			my_en_int_RRR	=	en_int_DGN_RRR( E1, E2, E4, smr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E2-E3)	<	kubo_tol ) 			&	!	E2==E3
				.and.	(	abs(E1-E2)	>	kubo_tol )			&	!	E1/=E2)then)then
	)then
		!	Eq.(B9):
			my_en_int_RRR	=	en_int_DGN_RRR(E2, E1, E4, smr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E2-E3)	>=	kubo_tol )			&	
				.and.	(	abs(E1-E3)	>=	kubo_tol )			&
				.and.	(	abs(E1-E2)	>=	kubo_tol )			&
	)then
		smrsmr	= smr**2
		!	EQ.(B7):	 I1(E1,E2,E3,E4)
		dE_34		=	E4 -E3 
		dE_24		=	E4 -E2 
		!
		denom_1213	=	(E1-E2) * (E1-E3)
		denom_2321	=	(E2-E3) * (E2-E1)
		denom_3132	=	(E3-E1) * (E3-E2)
		!
		!
		my_en_int_RRR	=	&
						log(	1.0_dp	+	dE_14**2	/smrsmr	)	/	(	2.0_dp	* denom_1213	)		&
					+	log(	1.0_dp	+	dE_24**2	/smrsmr	)	/	(	2.0_dp	* denom_2321	)		&
					+	log(	1.0_dp	+	dE_34**2	/smrsmr	)	/	(	2.0_dp	* denom_3132	)		&
					!
					+	atan(	dE_14			/  smr	)	/	(	i_dp 	* denom_1213	)				&
					+	atan(	dE_24			/  smr	)	/	(	i_dp 	* denom_2321	)				&
					+	atan(	dE_34			/  smr	)	/	(	i_dp 	* denom_3132	)											
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	else
		my_en_int_RRR	=	0.0_dp
	end if
	!
	my_en_int_RRR		=	my_en_int_RRR	
	!
	return
end function 
!~
!~
complex(dp) pure function en_int_DGN_RRR(E1, E3, E4,  smr)
	real(dp),		intent(in)		::	E1,E3,E4,smr
	real(dp)						::	dE_13, dE_14, dE_34
	!		see	 Freimuth et al., PRB 94, 144432 (2016)  
	!		 EQ.(B9)
	!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	dE_13		=	E3-E1
	dE_14		=	E4-E1
	dE_34		=	E4-E3
	!
	!
	en_int_DGN_RRR	=		log( 	( smr**2 + dE_34**2	)	/( smr**2 + dE_14**2 )	)			/	(2.0_dp * (dE_13)**2	)&
						- 	i_dp* 	(	atan(dE_34/smr)		-	atan(dE_14/smr)		)			/	(		dE_13**2		)&
						+	1.0_dp																/	(dE_13 * (dE_14+i_dP*smr))
	!
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
complex(dp) pure function my_en_int_RRA(E1, E2, E3, E4,  smr)
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
			my_en_int_RRA	=	&																														!
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
		my_en_int_RRA	=	&																															!
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
	!
	return
end function 
!~
!~


end module 