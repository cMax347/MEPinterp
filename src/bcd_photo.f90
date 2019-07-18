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
	public					::		bcd_photo_scnd
	!
contains



	subroutine bcd_photo_scnd(en_k, V_ka, hw_lst, smr_lst, ef_lst)
		real(dp),		intent(in)				::	hw_lst(:), smr_lst(:), ef_lst(:), en_k(:)
		complex(dp),	intent(in)				::	V_ka(:,:,:)
		real(dp),		allocatable				::	fd_k_deriv(:,:,:)
		real(dp)								::	fd_k_deriv_TOL=1e-3_dp
		!
		write(*,*)	"[bcd_photo_scnd]: WARNING fd_k_deriv_TOL hardcoded to ",fd_k_deriv_TOL
		!
		allocate(	fd_k_deriv(3,size(ef_lst,1), size(V_ka,3)	))
		!
		fd_k_deriv	=	get_fd_k_deriv(en_k, V_ka, ef_lst, fd_k_deriv_TOL)

		return
	end subroutine





	pure function kernel_jerk_term(k, i,j, w1, w2, V_ka, scnd_fermi_deriv, smr_lst) result(xhi_jerk)
		integer,		intent(in)		::	k, i,j 
		real(dp),		intent(in)		::	w1, w2
		complex(dp),	intent(in)		::	V_ka(:,:,:)
		real(dp),		intent(in)		::	scnd_fermi_deriv(:,:,:,:), smr_lst(:)
		complex(dp)						::	xhi_jerk
		!
		xhi_jerk	=	cmplx(0.0_dp, 0.0_dp, dp)
		!
		!do n=1 ,size(V_ka,3)
		!	!kernel_jerk_term	=	kernel_jerk_term	+	V_ka(k,n,n) * scnd_fermi_deriv(i,j,n,#ef_idx)
		!end do
		!
		return
	end function


	pure function kernel_bcd_term(fd_deriv)	result(xhi_bcd)
		real(dp),	intent(in)		::	 fd_deriv(:,:)			!	fd_deriv(1:N_ef,1:N_wf)
		complex(dp)					::	xhi_bcd
		

		!allocate(		tmp_ij12(3, n_smr) )
		!do n =1, size(V_ka,3)
		!	do m=1, size(V_ka,3)
		!		dE_nm		=	en_k(m)	-	en_k(n)
		!		!
		!		re_denom	=	
		!		!
		!		do smr =	1, size(smr_lst)	
		!			tmp(::,smr) =	tmp(::,smr)	+ 	V_ka(:,)
!
!		!		end do
!
!
!		!		!
!		!	end do
		!end do

		return
	end function


	pure function get_fd_k_deriv(en_k, V_ka, ef_lst, fd_k_deriv_TOL) result(fd_k_deriv)
		!	arXiv:1907.02532v1 -> eq.(43)
		real(dp),		intent(in)				::	fd_k_deriv_TOL, ef_lst(:), en_k(:)
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
					fd_k_deriv(:,ef_idx,n)	=	V_ka(:,n,n)
				end if
				!
			end do
		end do
		!
		!
		return
	end function	


end module bcd_photo






















