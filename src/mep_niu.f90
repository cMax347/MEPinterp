module mep_niu
	use constants,		only:		dp
	use input_paras,	only:		kubo_tol,	valence_bands
	use matrix_math,	only:		my_Levi_Civita,			& 
									convert_tens_to_vect

	implicit none

	private
	public			::		mep_niu_CS,		&
							mep_niu_IC,		&
							mep_niu_LC	

	save

contains

		!
!MEP RESPONSES
	pure function mep_niu_CS(A_ka, Om_kab) result(cs_tens)
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
			do n0 = 1, valence_bands
				call convert_tens_to_vect(	Om_kab(:,:,n0,n0),		om_vect(:)	)
				cs_scal	= cs_scal +	 0.5_dp  *	dreal(	 dot_product(	A_ka(:,n0,n0)	,	om_vect(:)	)			)
			end do
			!
			do i = 1, 3
				cs_tens(i,i)	= real(cs_scal,dp)
			end do
		end if
		!
		return
	end function


	pure function mep_niu_LC(velo, En) result(F3)
		!	LOCAL KUBO CONTRIBUTION
		!
		complex(dp),		intent(in)		::	velo(:,:,:)
		real(dp),			intent(in)		::	En(:)
		real(dp)							::	F3(3,3)
		integer								::	n0, n, neglected, tot, 		&
												j, k, l
		real(dp) 							::	pre_fact, velo_nom(3), en_denom
		!
		neglected	= 0
		tot 		= 0
		!
		F3	= 0.0_dp
		do n0 = 1, valence_bands
			!
			!MIXING
			do n = 1, size(velo,2)
			 	if( n/= n0 )	then	
			 		en_denom	=	(	en(n0) - en(n)	)**3
			 		tot			= 	tot + 1
			 		if(abs(en_denom) > kubo_tol) then 
			 			!
						!TRIPLE PRODUCT			 		
			 			do j = 1, 3
			 				do k = 1, 3
			 					do l = 1, 3
			 						pre_fact	=	real(			my_Levi_Civita(j,k,l)						,dp)	
			 						!
			 						!	
			 						if(		 abs(pre_fact) 	> 1e-1_dp	) then
			 							velo_nom(:)	=	real(		velo(:,n0,n) * velo(k,n,n0) * velo(l,n0,n0)		, dp)
			 							F3(:,j)		= 	F3(:,j) 	+	pre_fact   * velo_nom(:)	/	en_denom	
			 						end if
			 						!
			 					end do
			 				end do
			 			end do
			 		else
			 			neglected = neglected + 1
			 		end if
			 	end if
			 	!
			 	!
			end do
		end do
		!
		return
	end function


	pure function mep_niu_IC(velo, En) result(F2)
		!	ITINERANT KUBO CONTRIBUTION
		!
		complex(dp),		intent(in)		::	velo(:,:,:)
		real(dp),			intent(in)		::	En(:)
		real(dp)							::	F2(3,3)
		integer								::	n0, m, n, 		&
												j, k, l,		&
												neglected, tot
		real(dp) 							::	pre_fact, velo_nom(3), en_denom
		!
		neglected	= 0
		tot 		= 0
		!
		F2	= 0.0_dp
		do n0 = 1, valence_bands
			!
			!MIXING
			do m = 1, size(velo,2)
				do n = 1, size(velo,2)
				 	if( n/= n0 .and. m/=n0 )	then
				 		en_denom	=	(	en(n0) - en(n)	)**2		 * 		(	en(n0) - en(m)	)  
				 		tot 		= 	tot + 1
				 		if( abs(en_denom) > kubo_tol) then
				 			!
				 			!TRIPLE PRODUCT
				 			do j = 1, 3
				 				do k = 1, 3
				 					do l = 1, 3
				 						pre_fact	= 	- real(my_Levi_Civita(j,k,l),dp)
				 						!
				 						!
				 						if(		 abs(pre_fact) 	> 1e-1_dp	) then
											velo_nom(:)	=	real(		velo(:,n0,n) * velo(k,n,m) * velo(l,m,n0)	, dp)
				 							F2(:,j)		= 	F2(:,j)		+	pre_fact  * velo_nom(:)	/	en_denom	
				 						end if
				 						!
				 					end do
				 				end do
				 			end do
				 		else
				 			neglected	= neglected + 1
				 		end if
				 	end if
				 	!
				 	!
				end do
			end do
		end do
		!
		return
	end function



end module mep_niu






!!privat old:
!	subroutine	getF2(nZero,ki, Velo ,En, F2)
!		!
!		!	F^(2)_ij = + Re \sum_{n/=0,m/=0} \eps_{j,k,l} * (V^k_nm V^l_m0 V^i_mn) / ( (E0-En)**2 (E0-Em) )
!		!
!		integer,		intent(in)		:: nZero, ki
!		complex(dp),	intent(in)		:: Velo(:,:,:,:)  
!		real(dp),		intent(in)		:: En(:,:)			!
!		real(dp),		intent(out)		:: F2(3,3)
!		complex(dp)						:: Vtmp
!		real(dp)						:: eDiff, eDiff1, eDiff2
!		integer							:: i, j, k, l, n,m, nSize
!		!
!		nSize	=	size(Velo,3)
!		F2		=	0.0_dp
!		!loop bands
!		do n = 1, nSize
!			if( n/=nZero ) then
!				do m = 1, nSize
!					if(  m/=nZero) then
!						!
!						!
!						!ENERGIES
!						eDiff1		= 	( 	En(nZero,ki) - En(n,ki)		)**2 
!						eDiff2		=  	( 	En(nZero,ki) - En(m,ki)		)
!						eDiff		= 	eDiff1 * eDiff2
!						!degenerate energy warning
!						if( abs(eDiff) < machineP )  then
!							write(*,*)	"[addF2]: WARNING for k point = ",ki
!							write(*,'(a,i3,a,i3,a,i3,a,e14.6)') "[addF2]: WARNING degenerate bands n0=",nZero,"n=",n," m=",m," eDiff=",eDiff
!							write(*,'(a,i3,a,i3,a,e14.6)')	"[addF2]: ( E(",nZero,")-E(",n,") )**2=", eDiff1
!							write(*,'(a,i3,a,i3,a,e14.6)')	"[addF2]: ( E(",nZero,")-E(",m,") )   =", eDiff2
!							write(*,*)	"[addF2]: E(nZero=",nZero,")=",En(nZero,ki)
!							write(*,*)	"[addF2]: E(n=",n,")=",En(n,ki)
!							write(*,*)	"[addF2]: E(m=",m,")=",En(m,ki)
!						end if		
!						!
!						!loop matrix indices
!						do j = 1, 3
!							do i = 1, 3
!								!
!								!loop levi civita
!								do k = 1, 3
!									do l = 1, 3				
!										!VELOCITIES
!										Vtmp		= Velo(k,n,m,ki) * Velo(l,m,nZero,ki) * Velo(i,nZero,n,ki) 
!										!if( dimag(Vtmp) > 1e-10_dp ) write(*,*)	"[addF2]: none zero imag velo product: ",dimag(Vtmp),"; real part: ",dreal(Vtmp)
!										!MATRIX
!										F2(i,j) 	= F2(i,j) +   real(my_Levi_Civita(j,k,l),dp) *  dreal( Vtmp )  / eDiff	
!										!F2(i,j) 	= F2(i,j) +   my_Levi_Civita(j,k,l) *  dreal( Vtmp )  / eDiff	
!									end do
!								end do
!								!
!							end do
!						end do
!						!
!					end if
!				end do
!			end if
!		end do
!		!
!		!
!		return
!	end subroutine
!
!
!
!
!	subroutine	getF3(prefactF3, nZero,ki, Velo ,En, F3)	
!		!
!		!	F^(2)_ij = +- Re \sum_{n/=0} \eps_{j,k,l}  * (v^k_0 V^l_nZero V^i_0n) / ( (E0-En)**3  )
!		!
!		integer,		intent(in)		:: nZero, ki
!		complex(dp),	intent(in)		:: Velo(:,:,:,:)  
!		real(dp),		intent(in)		:: prefactF3, En(:,:)			
!		real(dp),		intent(out)		:: F3(3,3)
!		complex(dp)						:: Vtmp
!		real(dp)						:: eDiff
!		integer							:: i, j, k, l, n, nSize
!		!
!		nSize 	=	size(Velo,3)
!		F3		=	0.0_dp
!		!loop bands
!		do n = 1, nSize
!			if( n/=nZero ) then
!				!ENERGIES
!				eDiff		= ( 	En(nZero,ki) - En(n,ki)	 )**3 
!				!degenerate energy warning
!				if( abs(eDiff) < machineP ) write(*,*) "[addF3]: WARNING degenerate bands n0=",nZero,"n=",n," eDiff=",eDiff
!				!
!				!loop matrix indices
!				do j = 1, 3
!					do i = 1, 3
!						!
!						!loop levi civita
!						do k = 1, 3
!							do l = 1,3				
!								!VELOCITIES
!								Vtmp		= Velo(k,nZero,nZero,ki) * Velo(l,n,nZero,ki) * Velo(i,nZero,n,ki) 
!								!if( dimag(Vtmp) > 1e-10_dp ) write(*,*)	"[addF3]: none zero imag velo product: ",dimag(Vtmp),"; real part: ",dreal(Vtmp)
!								!
!								!MATRIX
!								F3(i,j) 	= F3(i,j) + real(prefactF3,dp) * real(my_Levi_Civita(j,k,l),dp) *	 dreal( Vtmp ) / eDiff
!								!F3(i,j) 	= F3(i,j) + prefactF3 * my_Levi_Civita(j,k,l) *	 dreal( Vtmp ) / eDiff
!							end do								!
!						end do
!						!
!					end do
!				end do
!				!
!			end if
!		end do
!		!
!		!
!		return
!	end subroutine
!
!
!
!	subroutine	getF3essin(prefactF3, nZero,ki, Velo ,En, F3)	
!		!
!		!	F^(2)_ij = +- Re \sum_{n/=0} \eps_{j,k,l}  * (v^k_0 V^l_nZero V^i_0n) / ( (E0-En)**3  )
!		!
!		integer,		intent(in)		:: nZero, ki
!		complex(dp),	intent(in)		:: Velo(:,:,:,:)  
!		real(dp),		intent(in)		:: prefactF3, En(:,:)			
!		real(dp),		intent(out)		:: F3(3,3)
!		complex(dp)						:: Vtmp
!		real(dp)						:: eDiff
!		integer							:: i, j, k, l, n, nSize
!		!
!		nSize 	=	size(Velo,3)
!		F3		=	0.0_dp
!		!loop bands
!		do n = 1, nSize
!			if( n/=nZero ) then
!				!ENERGIES
!				eDiff		= ( 	En(nZero,ki) - En(n,ki)	 )**3 
!				!degenerate energy warning
!				if( abs(eDiff) < machineP ) write(*,*) "[addF3]: WARNING degenerate bands n0=",nZero,"n=",n," eDiff=",eDiff
!				!
!				!loop matrix indices
!				do j = 1, 3
!					do i = 1, 3
!						!
!						!loop levi civita
!						do k = 1, 3
!							do l = 1,3				
!								!VELOCITIES
!								Vtmp		= Velo(i,nZero,n,ki)  * Velo(k,n,nZero,ki) * Velo(l,nZero,nZero,ki)
!								!if( dimag(Vtmp) > 1e-10_dp ) write(*,*)	"[addF3]: none zero imag velo product: ",dimag(Vtmp),"; real part: ",dreal(Vtmp)
!								!
!								!MATRIX
!								F3(i,j) 	= F3(i,j) + real(prefactF3,dp) * real(my_Levi_Civita(j,k,l),dp) *	 dreal( Vtmp ) / eDiff
!								!F3(i,j) 	= F3(i,j) + prefactF3 * my_Levi_Civita(j,k,l) *	 dreal( Vtmp ) / eDiff
!							end do								!
!						end do
!						!
!					end do
!				end do
!				!
!			end if
!		end do
!		!
!		!
!		return
!	end subroutine
