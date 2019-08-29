module keldysh
	!

	!
	use omp_lib
	use input_paras,	only:	kubo_tol, E_int_MIN, N_energy_int
	use constants,		only:	dp, i_dp, pi_dp, aUtoEv	
	use statistics,		only:	fd_stat
	use matrix_math,	only:	mat_tr
	use green_en_int,	only:	en_int_RRR, en_int_RRA

	implicit none


	private
	public			::		keldysh_scnd_photoC, keldysh_scnd_photoC_NUMERICAL

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
function keldysh_scnd_photoC(kpt,	en_k, V_ka, hw_lst, smr_lst, ef_lst)			result(phi)
	!		see	 Freimuth et al., PRB 94, 144432 (2016)  
	!		 EQ.(2)
	real(dp),		intent(in)		::	kpt(3), en_k(:), hw_lst(:), smr_lst(:), ef_lst(:)
	complex(dp),	intent(in)		::	V_ka(:,:,:)
	complex(dp),	allocatable		::	phi(:,:,:,:,:,:)
	complex(dp)						::	v_ijk(3,3,3), v_ikj(3,3,3)
	integer							::	a,b,c, j,k, hw, smr, ef
	real(dp)						::	magic_kcut
	!
	allocate(	phi(	3,3,3,	size(hw_lst), size(smr_lst), size(ef_lst)		))
	phi	= cmplx(0.0_dp,0.0_dp,dp)
	!
	!	TODO:	OMP REDUCTION OVER "phi"
	!if(maxval(en_k)> (maxval(ef_lst)+maxval(hw_lst)))&	
	!	write(*,*)	'[keldysh_scnd_photoC]: WARNING bands beyond max relevant energy kpt=',kpt

	!magic_kcut	=	

	if( 	abs(	minval(en_k)-maxval(ef_lst)-maxval(hw_lst))	< 1e-6_dp )then
		magic_kcut=	norm2(kpt)
		write(*,*)	'[keldysh_scnd_photoC]:	minval(en_k)=',minval(en_k)*aUtoEv," (eV) at kpt=",kpt," (magic_kcunt=",magic_kcut,")" 
	end if



	!$OMP PARALLEL DO REDUCTION(+:phi) 									&
	!$OMP DEFAULT(NONE)	SHARED(V_ka, en_k, ef_lst, smr_lst, hw_lst)		&
	!$OMP PRIVATE(c,b,a, k,j, v_ijk, v_ikj, ef, smr, hw)
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
	!$OMP END PARALLEL DO
	!
	return
end function
!
!




function keldysh_scnd_photoC_NUMERICAL(	en_k, U_k, V_ka, hw_lst, smr_lst, ef_lst)			result(phi)
	real(dp),		intent(in)		::	en_k(:), hw_lst(:), smr_lst(:), ef_lst(:)
	complex(dp),	intent(in)		::	U_k(:,:), V_ka(:,:,:)
	complex(dp),	allocatable		::	phi(:,:,:,:,:,:), G_base(:,:,:), tmp(:,:), G_R(:,:), G_R_shift(:,:),G_A(:,:)
	integer							::	n, ef, int_point, smr, hw, k, j, i 
	real(dp)						::	epsilon
	!
	allocate(	phi(	3,3,3,	size(hw_lst), size(smr_lst), size(ef_lst)		))
	phi	= cmplx(0.0_dp,0.0_dp,dp)	
	!
	allocate(	G_base(			size(U_k,1),size(U_k,2),	size(en_k,1)	))
	allocate(	G_R(			size(U_k,1),size(U_k,2)						))
	allocate(	G_A(			size(U_k,1),size(U_k,2)						))
	allocate(	G_R_shift(		size(U_k,1),size(U_k,2)						))	
	allocate(	tmp(			size(U_k,1),size(U_k,2)						))
	!
	!	get |kn><kn| Matrix	  
	do n = 1, size(U_k,1)
		G_base(:,:,n)	=	matmul(U_k(:,n:n),U_k(n:n,:))
	end do
	!
	!
	!
	!$OMP PARALLEL DO  REDUCTION(+:phi) COLLAPSE(2)	&
	!$OMP DEFAULT(NONE)	SHARED(ef_lst,smr_lst, hw_lst,  en_k, V_ka, G_base, N_energy_int, E_int_MIN)	&
	!$OMP PRIVATE(ef, int_point, epsilon, smr, G_R, G_A, hw, G_R_shift, k,j,i, tmp )
	do ef = 1, size(ef_lst)	
		do int_point = 1, N_energy_int
			epsilon	=	E_int_MIN		 + (ef_lst(ef)-E_int_MIN) 	*	 real(int_point-1,dp)/real(N_energy_int-1,dp)	
			!
			!
			!	TREADS I-IV (	all terms with f(epislon) prefactor )
	 		if(epsilon< ef_lst(ef))	then
	 			do smr = 1 , size(smr_lst)
	 				G_R			=	get_Green_retarded(			epsilon			, 	en_k, G_base, smr_lst(smr)	)
	 				G_A			=	conjg(transpose(G_R))
	 				!
	 				!	I & II		 contribution
					do hw =1, size(hw_lst)
	 					G_R_shift	=	get_Green_retarded(		epsilon-hw_lst(hw)	,	en_k, G_base, smr_lst(smr)	)
	 					do k =1,3
	 						do j=1,3
	 							do i=1,3
	 								tmp(:,:)	=	velo_propagator(			i,j,k,			V_ka, G_R, G_R_shift)
	 								!
	 								phi(i,j,k,hw,smr,ef)	=		phi(i,j,k,hw,smr,ef)	+	mat_tr(	matmul(	tmp, 	G_R	)	)
	 								phi(i,j,k,hw,smr,ef)	=		phi(i,j,k,hw,smr,ef)	-	mat_tr( matmul(	tmp, 	G_A	)	)
	 							end do
	 						end do
	 					end do
	 				end do
	 				!
					!	III & IV		 contribution
	 				do hw =1, size(hw_lst)
						G_R_shift	=	get_Green_retarded(		epsilon+hw_lst(hw),		en_k, G_base, smr_lst(smr)	)
						do k =1,3
	 						do j=1,3
	 							do i=1,3
	 								tmp(:,:)	=	velo_propagator(			i,k,j,			V_ka, G_R, G_R_shift)
	 								!
	 								phi(i,j,k,hw,smr,ef)	=		phi(i,j,k,hw,smr,ef)	+	mat_tr(	matmul(	tmp, 	G_R	)	)	
	 								phi(i,j,k,hw,smr,ef)	=		phi(i,j,k,hw,smr,ef)	-	mat_tr(	matmul(	tmp, 	G_A	)	)
	 							end do
	 						end do
	 					end do
	 				end do
	 			end do
	 		end if
	 		!
	 		!
			do smr = 1 , size(smr_lst)
	 			G_R			=	get_Green_retarded(			epsilon			, 	en_k, G_base, smr_lst(smr)	)
	 			G_A			=	conjg(transpose(G_R))
	 			!
	 			do hw =1, size(hw_lst)
	 				!
	 				!	TREADS V (	f(epislon-hw) term )
					if(	epsilon-hw_lst(hw) < ef_lst(ef)) then						
	 					G_R_shift	=	get_Green_retarded(		epsilon-hw_lst(hw),		en_k, G_base, smr_lst(smr)	)
	 					do k =1,3
	 						do j=1,3
	 							do i=1,3
									tmp(:,:)				=	velo_propagator(	i,j,k,			V_ka, G_R, G_R_shift)
									!
									phi(i,j,k,hw,smr,ef)	=	phi(i,j,k,hw,smr,ef)	+		mat_tr(matmul(	tmp, 	G_A	))
								end do	
							end do
						end do
	 				end if
	 				!
	 				!	TREADS VI (	f(epislon+hw) term )
	 				if(	epsilon+hw_lst(hw) < ef_lst(ef)) then
	 					G_R_shift	=	get_Green_retarded(		epsilon+hw_lst(hw),			 en_k, G_base, smr_lst(smr)	)
						do k =1,3
	 						do j=1,3
	 							do i=1,3
									tmp(:,:)	=	velo_propagator(			i,k,j,						V_ka, G_R, G_R_shift	)
									!
									phi(i,j,k,hw,smr,ef)	=		phi(i,j,k,hw,smr,ef)	+	mat_tr(matmul(	tmp, 	G_A	))
								end do
							end do
						end do	
	 				end if
	 			end do
	 		end do
	 		!
	 		!
	 	end do
	end do
	!$OMP END PARALLEL DO 
	!
	!	NORMALIZE ENERGY INTEGRATION (TODO: DO WE NEED MORE HERE)
	phi	=	phi	/	N_energy_int

	!
	!
end function


pure function velo_propagator(a,b,c,	V_ka, G_R, G_R_shift)	result(V_abc)
	integer,		intent(in)			::	a,b,c
	complex(dp),	intent(in)			::	V_ka(:,:,:), G_R(:,:), G_R_shift(:,:)
	complex(dp)	,	allocatable			::	V_abc(:,:)
	!
	allocate(	V_abc(	size(G_R,1),size(G_R,2)		))
	!
	V_abc(:,:)	=	matmul(	V_ka(a,:,:),	G_R(:,:)		)
	V_abc(:,:)	=	matmul(	V_abc(:,:),		V_ka(b,:,:)		)
	V_abc(:,:)	=	matmul(	V_abc(:,:),		G_r_shift(:,:)	)
	V_abc(:,:)	=	matmul(	V_abc(:,:),		V_ka(c,:,:)		)	
	!
	return
end function	






pure function get_Green_retarded(	epsilon, en_k, G_base, 	smr_val)		result(G_R)
	real(dp),		intent(in)		::	epsilon, en_k(:), smr_val
	complex(dp),	intent(in)		::	G_base(:,:,:)
	complex(dp),	allocatable		::	G_R(:,:)
	integer							::	n
	!
	allocate(	G_R(	size(G_base,1),size(G_base,2)		))
				G_R	=	cmplx(0.0_dp,0.0_dp,dp)
	!
	!
	do n = 1, size(en_k,1)
		G_R(:,:)	=		G_R(:,:)	+	G_base(:,:,n)	/	(epsilon - en_k(n) + i_dp*smr_val)
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














end module keldysh










