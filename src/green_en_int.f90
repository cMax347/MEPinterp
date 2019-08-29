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
	public 		::				en_int_RRR,					&
								en_int_RRA
	!
	real(dp),	parameter	::	pi_half	=	pi_dp / 2.0_dp, &
								err_tol	=	1e-1_dp
	
	logical,	parameter	::	verbose_RRR	=	.True.,	&
								verbose_RRA	=	.False.
	!
contains




	complex(dp)  function en_int_RRR(E1, E2, E3, E4, smr)
		real(dp),		intent(in)		::	E1, E2, E3, E4, smr
		complex(dp)							::	frank_RRR 
		complex(dp)						::	my_RRR
		!
		!	call my routine
		my_RRR	=	my_en_int_RRR(E1,E2,E3, E4, smr)	
		!
		!!	call franks routine
		!call integrals_rrr(E1,E2,E3, E4,smr,	frank_RRR)
		!call integrals_rrr(real(E1),real(E2),real(E3), real(E4),real(smr),	frank_RRR)
		!	compare and return franks
		
		!if(abs(my_RRR-frank_RRR)> err_tol	.and. verbose_RRR ) then
		!	!write(*,'(a,e10.3,a,e10.3,a,f6.2,a,f6.2,a,f6.2,a,f6.2)')					&	
		!			write(*,*)	'[en_int_RRR]	WARNING ',frank_RRR,' <-> ',my_RRR
		!end if 
		en_int_RRR 	=	my_RRR!
		!
		return
	end function


	complex(dp)  function en_int_RRA(E1, E2, E3, E4, smr)
		real(dp),		intent(in)		::	E1, E2, E3, E4, smr
		complex(dp)							::	frank_RRA 
		complex(dp)						::	my_RRA
		real(dp)						::	abs_err
		!
		!	call my routine
		my_RRA	=	my_en_int_RRA(E1,E2,E3, E4, smr)
		!	call franks routine
		!call integrals_rra(E1,E2,E3, E4, smr,	frank_RRA)
		!	compare and return franks
		!abs_err		=	abs(my_RRA-frank_RRA)

		!if( abs_err > err_tol .and. verbose_RRA ) then
		!	!write(*,*)	'[en_int_RRA]	WARNING ',frank_RRA,' <-> ',my_RRA
		!	write(*,*)	'[en_int_RRA]	WARNING (f_RRA=',frank_RRA,'<->',my_RRA,')'
		!else
		!!	write(*,'(a,e10.3)')		"[en_int_RRA]: 	accepted for dE_12=",E2-E1
		!!write(*,*)	'^^^^^^'
		!!write(*,*)	'[en_int_RRA]	sol=',frank_RRA,' at (',E1,', ',E2,', ', E3,')'
		!end if
		en_int_RRA 	=	my_RRA!
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
	real(dp)						::	denom_1213, denom_2321, denom_3132, dE_14, dE_24, dE_34,&
										re_rrr, im_rrr
	!
	!
	dE_14		=	E4 -E1
	smrsmr		= 	smr**2
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
		!	Eq.(B11):
		my_en_int_RRR	=	-0.5_dp	/	cmplx(dE_14, smr,dp)**2	
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E1-E2)	<	kubo_tol ) 			&	!	E1==E2
				.and.	(	abs(E1-E3)	>	kubo_tol )			&	!	E1/=E3
	)then
		!	Eq.(B9):
			my_en_int_RRR	=	en_int_DGN_RRR(	E1, E3, E4, smr, smrsmr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E1-E3)	<	kubo_tol ) 			&	!	E1==E3
				.and.	(	abs(E2-E1)	>	kubo_tol )			&	!	E1/=E2)then
	)then
		!	Eq.(B9):
			my_en_int_RRR	=	en_int_DGN_RRR( E1, E2, E4, smr, smrsmr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E2-E3)	<	kubo_tol ) 			&	!	E2==E3
				.and.	(	abs(E1-E2)	>	kubo_tol )			&	!	E1/=E2)then)then
	)then
		!	Eq.(B9):
			my_en_int_RRR	=	en_int_DGN_RRR(E2, E1, E4, smr, smrsmr)
		!________________________
	!^^^^^^^^^^^^^^^^^^^^^^^^^
	else if(			(	abs(E2-E3)	>=	kubo_tol )			&	
				.and.	(	abs(E1-E3)	>=	kubo_tol )			&
				.and.	(	abs(E1-E2)	>=	kubo_tol )			&
	)then
		!	EQ.(B7):	 I1(E1,E2,E3,E4)
		dE_34		=	E4 -E3 
		dE_24		=	E4 -E2 
		!
		denom_1213	=	(E1-E2) * (E1-E3)
		denom_2321	=	(E2-E3) * (E2-E1)
		denom_3132	=	(E3-E1) * (E3-E2)
		!
		!
		re_rrr		=	0.5_dp * log(	1.0_dp	+	dE_14**2	/smrsmr	)	/	denom_1213			&
					  +	0.5_dp * log(	1.0_dp	+	dE_24**2	/smrsmr	)	/	denom_2321			&
					  +	0.5_dp * log(	1.0_dp	+	dE_34**2	/smrsmr	)	/	denom_3132			
		!
		im_rrr		=	- 	atan(	dE_14			/  smr	)				/	denom_1213			&
						-	atan(	dE_24			/  smr	)				/	denom_2321			&
						-	atan(	dE_34			/  smr	)				/	denom_3132							
		!
		!
		my_en_int_RRR	=	cmplx(	re_rrr,	im_rrr,		dp)			
		!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	else
		my_en_int_RRR	=	cmplx(0.0_dp,0.0_dp,dp)
	end if
	!
	return
end function 
!~
!~
complex(dp) pure function en_int_DGN_RRR(E1, E3, E4,  smr, smrsmr)
	real(dp),		intent(in)		::	E1,E3,E4,smr, smrsmr
	real(dp)						::	dE_13, dE_14, dE_34, dE_1313
	!		see	 Freimuth et al., PRB 94, 144432 (2016)  
	!		 EQ.(B9)
	!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	dE_13		=	E3-E1
	dE_1313		=	dE_13**2
	dE_14		=	E4-E1
	dE_34		=	E4-E3
	!
	!
	en_int_DGN_RRR	=	cmplx(		0.5_dp * log( 	(  smrsmr + dE_34**2 )/( smrsmr + dE_14**2 )		)	/	 dE_1313	,&
									- 	 	(	atan(dE_34/smr)		-	atan(dE_14/smr)		)				/	dE_1313		,&
						dp)
	!
	en_int_DGN_RRR	=	en_int_DGN_RRR	+	cmplx(1.0_dp,0.0_dp,dp)	 /	( cmplx(dE_13,0.0_dp,dp) * cmplx(dE_14,smr,dp)	)
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
complex(dp)  function my_en_int_RRA(E1, E2, E3, E4,  smr)
	!		see	 Freimuth et al., PRB 94, 144432 (2016) 
	!		 EQ.(B8)	/	EQ.(B10)	
	!		i.e.	I2(E1,E2,E3,E4)
	real(dp),		intent(in)		::	E1, E2, E3, E4, smr
	real(dp)						::	dE_12,	dE_13,	dE_14,	&
											  	dE_23,	dE_24,	&
											  			dE_34,	&
										rval14, rval24, rval34	  			

	complex(dp)						::	cnom14, cnom24, cnom34,		&
										cdenom_13, cdenom_1313,  	&
										cdenom_1213, cdenom_2312, cdenom_1323
	!
	!	
	dE_13		=	E3-E1
	dE_14		=	E4-E1
	dE_34		=	E4-E3
	!
	cdenom_13 		=	cmplx(dE_13,2.0_dp*smr,dp)		
	cdenom_1313 	=	cdenom_13**2
	!
	if(		abs(E1-E2)	<	kubo_tol )	then
		!
		!write(*,*)	'[my_en_int_RRA]: E1 E2 degenazryc'
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
		!	Eq.(B10):			I2(E1,E1,E3,E4)																										!
			my_en_int_RRA	=	cmplx(0.5_dp * log( (smr**2+dE_34**2)/	(smr**2+dE_14**2)), 0.0_dp,dp)		/		cdenom_1313		&							!
						 	 +	cmplx(0.0_dp,		pi_half +atan(dE_34	/ smr)	,dp)						/		cdenom_1313		&							!
						 	 +	cmplx(0.0_dp,		pi_half	+atan(dE_14 / smr)	,dp)						/		cdenom_1313									!
			!~
			my_en_int_RRA	= my_en_int_RRA		 	 +	 cmplx(1.0_dp,0.0_dp,dp)	/	( cdenom_13*cmplx(dE_14,smr,dp)	)							!
		!___________________________________________________________________________________________________________________________________________!
	else
		!
		dE_12		=	E2-E1
		dE_23 		=	E3-E2
		dE_24		=	E4-E2
		!
		rval14			=	1.0_dp+(dE_14/smr)**2
		rval24			=	1.0_dp+(dE_24/smr)**2
		rval34			=	1.0_dp+(dE_34/smr)**2
		!
		!write(*,*)	rval14,' ',rval24,' ',rval34
		!	setup the nominators
		cnom14		=	cmplx(		0.5_dp*log(rval14),			pi_half+ atan(dE_14/smr) 		,dp)
		cnom24		= -	cmplx( 		0.5_dp*log(rval24),			pi_half+ atan(dE_24/smr)		,dp)	
		cnom34		= 	cmplx(		0.5_dp*log(rval34),			pi_half+ atan(dE_34/smr)		,dp)
		!
		!	setup denominators
		cdenom_1213	=		cmplx(	dE_12		,      0.0_dp	,dp) 	* cdenom_13	
		cdenom_2312	=		cmplx(dE_23*dE_12	, 	2.0_dp*smr	,dp)
		cdenom_1323	=		cmplx(	dE_23		, 	2.0_dp*smr	,dp)	* cdenom_13 
		!
		!
		!	build the sum
		my_en_int_RRA	=						cnom14	/	cdenom_1213
		my_en_int_RRA	=	my_en_int_RRA	+	cnom24	/	cdenom_2312
		my_en_int_RRA	=	my_en_int_RRA	+	cnom34	/	cdenom_1323

		!write(*,*)	'[my_en_int_RRA]:  tot=	',my_en_int_RRA	!/	cdenom_1213
		!write(*,*)	'[my_en_int_RRA]: sum1=	',cnom14	!/	cdenom_1213
		!write(*,*)	'[my_en_int_RRA]: sum2=	',cnom24	!/	cdenom_2312
		!write(*,*)	'[my_en_int_RRA]: sum3=	',cnom34	!/	cdenom_1323
		!write(*,*)	'[my_en_int_RRA]:dE_12=	',dE_12
		!write(*,*)	'~'
		!___________________________________________________________________________________________________________________________________________!
	end if
	!
	return
end function 
!~
!~
!complex(dp) pure function my_en_int_RRA(E1, E2, E3, E4,  smr)
!	!		see	 Freimuth et al., PRB 94, 144432 (2016) 
!	!		 EQ.(B8)	/	EQ.(B10)	
!	!		i.e.	I2(E1,E2,E3,E4)
!	real(dp),		intent(in)		::	E1, E2, E3, E4, smr
!	complex(dp)						::	i_2smr,smrsmr
!	!
!	!	
!	i_2smr		=	2.0_dp * i_dp * smr	
!	smrsmr		=	smr**2	
!	!
!	if(		(abs(E1-E2)	< kubo_tol)		)then
!		!
!		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!		!	Eq.(B10):			I2(E1,E1,E3,E4)																										!
!			my_en_int_RRA	=	&																													!
!						 + log( 	(smrsmr+(E3-E4)**2)/	(smrsmr+(E1-E4)**2))/		(	2.0_dp* (E1-E3-i_2smr)**2	)	&						!
!						 + i_dp * (		pi_half 	+	atan((E4-E3)/smr)	)	/		(			(E1-E3-i_2smr)**2	)	&						!
!						 + i_dp * (		pi_half		+	atan((E4-E1)/smr)	)	/		(			(E1-E3-i_2smr)**2	)	&						!
!						 +  1.0_dp												/		((E3-E1+i_2smr)*(E4-E1+i_dp*smr))							!
!		!___________________________________________________________________________________________________________________________________________!
!	else
!		!
!		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!		!	Eq.(B8):				I2(E1,E2,E3,E4)																									!
!																																					!																																	!
!		my_en_int_RRA	=	&																														!
!					+	log(	1.0_dp	+	(E1 - E4)**2 	/smrsmr		)		/	(	2.0_dp *	 (E1-E2) 		*	( E1-E3 -i_2smr)		)&	!		
!					+	log(	1.0_dp	+	(E2 - E4)**2	/smrsmr		)		/	(	2.0_dp * ( E2-E3 -i_2smr) 	* 		(E2-E1)				)&	!
!					+	log(	1.0_dp	+	(E3 - E4)**2	/smrsmr		)		/	(	2.0_dp * ( E3-E1 +i_2smr) 	*	( E3-E2 +i_2smr)		)&	!
!					!-----																															!
!					+	i_dp * (	pi_half		+	atan((E4-E1)/smr)	)		/	(			(E1-E2) 			*	( E3-E1 + i_2smr)		)&	!
!					+	i_dp * (	pi_half		+	atan((E4-E2)/smr)	)		/	(		( E3-E2 + i_2smr)		*		(E2-E1)				)&	!
!					+	i_dp * (	pi_half		+	atan((E4-E3)/smr)	)		/	(		( E3-E1 + i_2smr)		*	( E3-E2 + i_2smr)		)	!					
!		!___________________________________________________________________________________________________________________________________________!
!	end if
!	!
!	!
!	return
!end function 


end module 