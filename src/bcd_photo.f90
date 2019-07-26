!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!				arXiv:1907.02532v1 
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module bcd_photo
	use constants,		only:		dp, i_dp, pi_dp
	use statistics,		only:		fd_stat
	use omp_lib
	!
	implicit none
	!
	private	
	public					::		photoC_bcd_rect_term,	&
									photoC_bcd_term,		& 
									photoC_jrk_term,		&
									photoC_inj_term,		&
									photoC_ib_term,			&
									get_fd_k_deriv
	!
	real(dp),	parameter	::		fd_k_deriv_TOL=1e-2_dp
	logical,	save		::		printed_jerk=.False.
	
	!
contains


!~~
!~~
!~~
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!		METALLIC CONTRIBUTIONS
!	
!~~~
	pure function photoC_bcd_term(hw_lst, smr_lst, ef_lst, fd_k_deriv,en_k, V_ka)	result(sigma_bcd_raw)
		real(dp),				intent(in)			::	hw_lst(:), smr_lst(:), ef_lst(:), fd_k_deriv(:,:,:), en_k(:)
		complex(dp),			intent(in)			::	V_ka(:,:,:)
		complex(dp),			allocatable			::	sigma_bcd_raw(:,:,:,:,:,:)
		complex(dp)									::	tmp_ij(3), tmp_ji(3)
		integer										::	ef, smr, hw, i, j
		!
		allocate(	sigma_bcd_raw( 3,3,3,		size(hw_lst), size(smr_lst), size(ef_lst)		))
		sigma_bcd_raw	=	cmplx(0.0_dp,0.0_dp,dp)
		!
		do ef=1, size(ef_lst)
			do smr=1, size(smr_lst)
				do hw =1, size(hw_lst)
					do i = 1,3
						do j = 1,3 
							!
							tmp_ij(:)	=	bcd_kernel(	-hw_lst(hw),	hw_lst(hw),		i,j, 		smr_lst(smr),ef, fd_k_deriv, en_k, V_ka)
							tmp_ji(:)	=	bcd_kernel(	 hw_lst(hw),  - hw_lst(hw),		j,i,		smr_lst(smr),ef, fd_k_deriv, en_k, V_ka)
							!
							sigma_bcd_raw(:,j,i,hw,smr,ef)	=		tmp_ij(:) + tmp_ji(:)
						end do
					end do 
				end do 
			end do 
		end do 
		!
		return
	end function
	!~
	!~
	function photoC_jrk_term(hw_lst, smr_lst, ef_lst, fd_k_deriv,en_k, V_ka)	result(sigma_jrk)
		real(dp),			intent(in)				::	hw_lst(:), smr_lst(:), ef_lst(:), fd_k_deriv(:,:,:), en_k(:) !	fd_deriv(1:N_ef,1:N_wf)
		complex(dp),		intent(in)				::	V_ka(:,:,:)
		complex(dp), 		allocatable				::	sigma_jrk(:,:,:,:,:,:)
		!
		allocate(	sigma_jrk(3,3,3,  size(hw_lst),	size(smr_lst), size(ef_lst))		)
		sigma_jrk	=	cmplx(0.0_dp, 0.0_dp, dp)
		!
		if(.not. printed_jerk) then
			write(*,*)	"[photoC_jrk_term]: WARNING not implemented yet"		
			printed_jerk	=	.True.
		end if
		!
		return
	end function
	!~
	!~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~
!~~
!~~
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!		INSULATING CONTRIBUTIONS
!	
!~~~
	pure function	photoC_INJ_term(hw_lst, smr_lst, fd_distrib, en_k, V_ka)	result(sigma_inj)
		real(dp),			intent(in)				::	hw_lst(:), smr_lst(:), fd_distrib(:,:), en_k(:) !	fd_deriv(1:N_ef,1:N_wf)
		complex(dp),		intent(in)				::	V_ka(:,:,:)
		complex(dp), 		allocatable				::	sigma_inj(:,:,:,:,:,:)
		complex(dp)									::	tmp_ij(3), tmp_ji(3)
		integer										::	ef, smr, hw, i, j
		!
		allocate(	sigma_inj(3,3,3,  size(hw_lst),	size(smr_lst), size(fd_distrib,1))	)
		sigma_inj	=	cmplx(0.0_dp,0.0_dp,dp)
		!
		do ef=1, size(fd_distrib,1)
			do smr=1, size(smr_lst)
				do hw =1, size(hw_lst)
					do i = 1,3
						do j = 1,3 
							!
							tmp_ij(:)	=	inj_kernel(	-hw_lst(hw),	hw_lst(hw),		i,j, 		smr_lst(smr),ef, fd_distrib, en_k, V_ka)
							tmp_ji(:)	=	inj_kernel(	 hw_lst(hw),  - hw_lst(hw),		j,i,		smr_lst(smr),ef, fd_distrib, en_k, V_ka)
							!
							sigma_inj(:,j,i,hw,smr,ef)	=		tmp_ij(:) + tmp_ji(:)
						end do
					end do 
				end do 
			end do 
		end do 
		!
		return
	end function
	!~	
	!~
	pure function photoC_ib_term(hw_lst, smr_lst, ef_lst, fd_k_deriv,en_k, V_ka, Om_kab)	result(sigma_ib)
		real(dp),			intent(in)				::	hw_lst(:), smr_lst(:), ef_lst(:), fd_k_deriv(:,:,:), en_k(:) !	fd_deriv(1:N_ef,1:N_wf)
		complex(dp),		intent(in)				::	V_ka(:,:,:), Om_kab(:,:,:,:)
		complex(dp), 		allocatable				::	sigma_ib(:,:,:,:,:,:)
		allocate(	sigma_ib(3,3,3,  size(hw_lst),	size(smr_lst), size(ef_lst))	)
		sigma_ib	=	cmplx(0.0_dp,0.0_dp,dp)
		!
		return
	end function
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








!~~
!~~
!~~
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!		PUBLIC HELPER
!	
!~~~
	pure function get_fd_k_deriv(en_k, V_ka, ef_lst) result(fd_k_deriv)
		!	arXiv:1907.02532v1 -> eq.(43)
		real(dp),		intent(in)				::	ef_lst(:), en_k(:)
		complex(dp),	intent(in)				::	V_ka(:,:,:)
		real(dp),		allocatable				::	fd_k_deriv(:,:,:)
		integer									::	n, ef_idx
		!
		allocate(	fd_k_deriv(3,size(ef_lst),size(V_ka,3))		) 
		fd_k_deriv	=	0.0_dp
		!
		do n=1, size(en_k,1)
			do ef_idx=1, size(ef_lst,1)
				!
				if(			abs(en_k(n) - ef_lst(ef_idx))	< fd_k_deriv_TOL 			) then
					!
					fd_k_deriv(:,ef_idx,n)	=	real(	V_ka(:,n,n), dp)					!
				end if
			end do
		end do
		!~
		!~
		return
	end function
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






!~~
!~~
!~~
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!		KERNELS
!	
!~~~
	pure function bcd_kernel(w1,w2,i,j, smr_val, ef_idx, fd_k_deriv, en_k, V_ka) result(sigma_bcd_vec)
		!^^^^^^^^		
		!
		!	see: Matsyshyn & Sodemann arXiv:1901.04467v1 eq.(27)
		!		https://arxiv.org/pdf/1907.02532.pdf
		!
		!^^^^^^^^
		real(dp),			intent(in)		::	w1,w2, en_k(:), smr_val,fd_k_deriv(:,:,:)
		integer,			intent(in)		::	i,j, ef_idx
		complex(dp),		intent(in)		::	V_ka(:,:,:)
		complex(dp)							::	sigma_bcd_vec(3)
		integer								::	m,n
		!
		sigma_bcd_vec	=	cmplx(0.0_dp,0.0_dp,dp)
		!
		do m = 1, size(en_k)
			do n = 1, size(en_k)
				if(m==n)cycle
				!
				sigma_bcd_vec(:)	=	sigma_bcd_vec(:)															&
						+	V_ka(:,m,n) * V_ka(i,n,m) * (	fd_k_deriv(j,ef_idx,m) - fd_k_deriv(j,ef_idx,n)		)	&	
						/ 	(	(w2+i_dp*smr_val) 	* 	(w1+w2+en_k(m)-en_k(n)+i_dp*smr_val)	 *	 (en_k(n)-en_k(m))		)
			end do
		end do
		! 
		sigma_bcd_vec	=	-0.5_dp * sigma_bcd_vec
		!
		return
	end function
	!~
	!~
	pure function inj_kernel(w1,w2,i,j, smr_val, ef_idx, fd_distrib, en_k, V_ka) result(sigma_inj_vec)
		!^^^^^^^^		
		!
		!	see: Matsyshyn & Sodemann arXiv:1901.04467v1 eq.(30)
		!		https://arxiv.org/pdf/1907.02532.pdf
		!
		!^^^^^^^^
		real(dp),			intent(in)		::	w1,w2, en_k(:), smr_val,fd_distrib(:,:)	!fd_distrib(1:N_ef,1:N_wf)
		integer,			intent(in)		::	i,j, ef_idx
		complex(dp),		intent(in)		::	V_ka(:,:,:)
		complex(dp)							::	sigma_inj_vec(3), freq_fact, pre_fact
		integer								::	m,n
		real(dp)							::	fd_fact
		!
		sigma_inj_vec	=	cmplx(0.0_dp,0.0_dp,dp)
		!
		!	SUM OVER STATES
		do m = 1 , size(en_k)
			do n = 1, size(en_k)
				if(n==m) cycle
				fd_fact		=	( fd_distrib(ef_idx,n) 	- 	fd_distrib(ef_idx,m) )	/	((	en_k(n)			-	en_k(m) )**2)		
				!
				if(abs(en_k(n)-en_k(m))>1e-15_dp) then
					freq_fact	=	fd_fact	/	(	(w1+en_k(n)-en_k(m)+i_dp*smr_val)	&
											  *		(w2-en_k(n)+en_k(m)+i_dp*smr_val)	)
					!
					sigma_inj_vec(:)	=	sigma_inj_vec(:)	+	freq_fact				& 	
												*(	(V_ka(:,n,n) 	- 	V_ka(:,m,m)	) 		&
													* V_ka(j,n,m)	*	V_ka(i,m,n)) 		
				end if
				!
			end do 
		end do
		!
		sigma_inj_vec	=	sigma_inj_vec * 0.5_dp * 	(	1.0_dp	 + 		i_dp* smr_val 			&
																		/ ( w1+w2+i_dp*smr_val)		&	
														)
		!	
		return
	end function
	!~
	!~
	pure function ib_kernel(w1,w2,i,j,	smr_val, ef_idx, fd_distrib,fd_k_deriv, en_k, V_ka,Om_kab)	result(sigma_ib_vec)
		!^^^^^^^^		
		!
		!	see: Matsyshyn & Sodemann arXiv:1901.04467v1 eq.(31)
		!		https://arxiv.org/pdf/1907.02532.pdf
		!
		!^^^^^^^^
		real(dp),			intent(in)		::	w1,w2, en_k(:), smr_val,fd_distrib(:,:), fd_k_deriv(:,:,:)	!fd_distrib(1:N_ef,1:N_wf)
		integer,			intent(in)		::	i,j, ef_idx
		complex(dp),		intent(in)		::	V_ka(:,:,:), Om_kab(:,:,:,:)
		complex(dp)							::	sigma_ib_vec(3), freq_fact, pre_fact, c_multi, c_sum(3), k_deriv
		integer								::	m,n, c
		real(dp)							::	fd_fact
		!
		sigma_ib_vec	=	cmplx(0.0_dp,0.0_dp,dp)
		!
		do m =1, size(en_k)
			do n=1, size(en_k)
				if(n==m)cycle
				!	add first term
				k_deriv	=	V_ka(j,n,m) 				* (fd_k_deriv(i,ef_idx,n)	-	fd_k_deriv(i,ef_idx,m)	)	/(en_k(m)-en_k(n))
				k_deriv = 	k_deriv + Om_kab(i,j,n,m)	* (fd_distrib(ef_idx,n)		-	fd_distrib(ef_idx,m)	)	
				!
				!
				sigma_ib_vec(:)		=	sigma_ib_vec(:) 	- i_dp * V_ka(:,m,n) * k_deriv 									&
											/ (	(w1+w2-en_k(n)+en_k(m)+i_dp*smr_val) * (w2-en_k(n)+en_k(m)+i_dp*smr_val)	)
				!~				
				!~
				!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
				!	sum over c term													
				c_sum		=	cmplx(0.0_dp,0.0_dp,dp)
				do c=1, size(en_k)
					c_sum(:)		=	c_sum(:) + 			V_ka(i,m,c) * V_ka(:,c,n) 									&
													/(	(en_k(c)-en_k(m)) * (w1+w2-en_k(n)+en_k(c)+i_dp*smr_val)		)
						!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
					c_sum(:)		=	c_sum(:) + 			V_ka(:,m,c) * V_ka(i,c,n) 									&
													/(	(en_k(n)-en_k(c)) * (w1+w2-en_k(c)+en_k(m)+i_dp*smr_val)		)	
						!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~							 
				end do
				c_multi				=	- ( fd_distrib(ef_idx,n)-fd_distrib(ef_idx,m) ) * V_ka(j,n,m)			&
											/ (		(en_k(m)-en_k(n))	*	(w2-en_k(n)+en_k(m)+i_dp*smr_val)	)
				c_sum(:)			=	c_sum(:) * c_multi
				!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				!
				!
				sigma_ib_vec(:)		=	sigma_ib_vec(:)	+ c_sum(:)
			end do
		end do
		!
		!	prefactor
		sigma_ib_vec		=	0.5_dp	* sigma_ib_vec
		!~
		!~		
		return
	end function








!~~
!~~
!~~
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!		KERNELS
!	
!~~~
	pure function photoC_bcd_rect_term(hw_lst, smr_lst, ef_lst, fd_k_deriv,en_k, V_ka)	result(sigma_bcd)
		!^^^^^^^^		
		!			SPECIAL CASE OF RECTIFICATION (w1=-w=-w2)
		!			
		!	THE AUTHOR OF THIS CODE DOES NOT TRUST EQ.(42) IMPLEMENTED IN THIS ROUTINE		!!!1!!!11!!!
		!	->	SINCE A DERIVATION FROM Eq.(27) (which is also used in above function bcd_kernel)
		!				
		!
		!
		!	see: Matsyshyn & Sodemann arXiv:1901.04467v1 eq.(42)
		!		https://arxiv.org/pdf/1907.02532.pdf
		!
		!^^^^^^^^
		real(dp),			intent(in)				::	hw_lst(:), smr_lst(:), ef_lst(:), fd_k_deriv(:,:,:), en_k(:) !	fd_deriv(1:N_ef,1:N_wf)
		complex(dp),		intent(in)				::	V_ka(:,:,:)
		complex(dp), 		allocatable				::	sigma_bcd(:,:,:,:,:,:)
		complex(dp)									::	velo_tens(3,3,3), &
														frq_broad, fermi_velo_j, cmplx_denom
		real(dp)									::	dE_nm, re_denom, im_denom, im_denom_smr, pre_fact
		integer										::	n, m, ef, smr, hw, k, j  
		!
		allocate(	sigma_bcd(3,3,3,  size(hw_lst),	size(smr_lst), size(ef_lst))		)
		sigma_bcd	=	cmplx(0.0_dp,0.0_dp,dp)
		pre_fact	=	- 0.5_dp
		!
		do n =1, size(V_ka,3)
			do m=1, size(V_ka,3)
				if(n==m)cycle
				!^^^^^
				dE_nm		=	en_k(n)	- en_k(m)
				!
				do ef = 1, size(ef_lst)
					!	APPLY SMEARING
					do smr = 1, size(smr_lst)
						cmplx_denom	=	1.0_dp / (	dE_nm *(-dE_nm+i_dp*smr_lst(smr))	)
						!~
						do hw = 1, size(hw_lst)
							frq_broad	=	pre_fact 	/	cmplx(hw_lst(hw),smr_lst(smr),dp)	
							!
							do k = 1, 3
								do j = 1,3 
									fermi_velo_j 		=	fd_k_deriv(j,ef,m) - fd_k_deriv(j,ef,n)	
									!
									sigma_bcd(:,j,k,hw,smr,ef)	 =		sigma_bcd(:,j,k,hw,smr,ef)					&
																	+ 	frq_broad * V_ka(:,n,m) * fermi_velo_j  * 	V_ka(k,m,n)   / cmplx_denom	
								end do
							end do
						end do	
					end do
				end do
				!~~
			end do
		end do
		!
		return
	end function









































end module bcd_photo






















