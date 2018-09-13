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
			!if( abs(imag(cs_scal))>1e-3_dp ) stop	'[mep_niu_CS]: ERROR found imaginary contributions'
			!	the CS is diagonal
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
												i, j, k, l
		complex(dp)							::	velo_nom
		real(dp) 							::	en_denom
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
			 				do i = 1, 3
			 					do l = 1, 3
			 						do k = 1, 3
			 							velo_nom	=	velo(i,n0,n) * velo(k,n,n0) * velo(l,n0,n0)
			 							!		
			 							F3(i,j)	= F3(i,j) + real(my_Levi_Civita(j,k,l),dp)  * real(	velo_nom	,dp )	/	en_denom	
			 							
			 						end do
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
		!if(neglected > 0) then
		!	write(*,'(a,i3,a,i5)',advance="no")		'[#',mpi_id,':get_F3]: dropped ',neglected
		!	write(*,'(a,i6,a)')						' of ',tot,' contributions due to degenerate bands'
		!end if
		return
	end function


	pure function mep_niu_IC(velo, En) result(F2)
		!	ITINERANT KUBO CONTRIBUTION
		!
		complex(dp),		intent(in)		::	velo(:,:,:)
		real(dp),			intent(in)		::	En(:)
		real(dp)							::	F2(3,3)
		integer								::	n0, m, n, 		&
												i, j, k, l,		&
												neglected, tot
		complex(dp)							::	velo_nom
		real(dp) 							::	en_denom
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
				 				do i = 1, 3
				 					do l = 1, 3
				 						do k = 1, 3
				 							velo_nom	=	velo(i,n0,n) * velo(k,n,m) * velo(l,m,n0)
				 							!
				 							F2(i,j)	= F2(i,j) - real(my_Levi_Civita(j,k,l),dp) * real(	velo_nom 	,dp )	/	en_denom	
				 						end do
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
		!if(neglected > 0) then
		!	write(*,'(a,i3,a,i5)',advance="no")		'[#',mpi_id,':get_F2]: dropped ',neglected
		!	write(*,'(a,i6,a)')						' of ',tot,' contributions due to degenerate bands'
		!end if
		return
	end function



end module mep_niu