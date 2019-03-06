module statistics
	use constants, 		only:			dp, aUtoEv, kBoltz_Eh_K
	!
	implicit none

	public					::			fd_stat,		fd_stat_deriv, 			&
										fd_get_N_el
	private


	interface fd_stat
		module procedure	d0_fd_stat
		module procedure	d1_fd_stat
	end interface fd_stat

	interface fd_stat_deriv 
		module procedure	d0_fd_stat_deriv
		module procedure	d1_fd_stat_deriv
	end interface fd_stat_deriv

	interface fd_get_N_el
		module procedure	d0_fd_get_N_el
		module procedure	d1_fd_get_N_el
	end interface fd_get_N_el		




	save


	real(dp),	parameter	::			min_temp	= 	1e-2_dp
	logical					::			t_zero_warn	=	.False.

contains





	real(dp) pure function d0_fd_stat(e_band, e_fermi,	T_kelvin) 
		real(dp), 		intent(in)			::	e_band, e_fermi, T_kelvin
		real(dp)							::	T_smear
		!
		d0_fd_stat		=	0.0_dp
		!
		if(	 T_kelvin > min_temp ) 				then
			!
			!	FINITE TEMPERATURE
			T_smear			=	kBoltz_Eh_K		*	T_kelvin
			d0_fd_stat		 	= 	1.0_dp	/	(	1.0_dp	+	exp(	(e_band	- e_fermi)	/	(T_smear)))
			!
			!
		else if(	e_band < e_fermi	)  		then
			!	
			!	ZERO TEMPERATURE
			d0_fd_stat	=	1.0_dp
		end if
		!
		!
		return
	end function


	pure function d1_fd_stat(e_band, ef_lst, T_kelvin)	result(fd_lst)
		real(dp),					intent(in)		::	e_band, ef_lst(:), T_kelvin
		real(dp),	allocatable						::	fd_lst(:)
		real(dp)									::	T_smear
		integer										::	ef_idx
		!
		allocate(	fd_lst(size(ef_lst,1))		)
		fd_lst	=	0.0_dp
		!
		if(T_kelvin > min_temp )	then
			!
			!	T/=0		FERMI-DIRAC SMEARING
			T_smear			=	kBoltz_Eh_K		*	T_kelvin
			do ef_idx	=	1, size(ef_lst,1)
				fd_lst(ef_idx)	=	1.0_dp /	(	1.0_dp + exp(	( e_band	- ef_lst(ef_idx) )	/	T_smear	))
			end do
		else
			!
			!	T==0		STEPFUNCTION
			do ef_idx	=	1, size(ef_lst,1)
				if(	e_band < ef_lst(ef_idx)	)		fd_lst(ef_idx)	=	1.0_dp
			end do
		end if
		!
		return
	end function





	pure function d0_fd_get_N_el(en_k, e_fermi, T_kelvin)	result(N_el)
		real(dp),		intent(in)			::	en_k(:), e_fermi, T_kelvin
		real(dp)							::	N_el
		integer								::	n
		!
		N_el	=	0.0_dp
		!
		do n = 1, size(en_k)
			N_el	=	N_el	+	fd_stat(en_k(n),	e_fermi, T_kelvin)
		end do
		!
		return
	end function


	pure function d1_fd_get_N_el(en_k, ef_lst, T_kelvin)	result(N_el_lst)
		real(dp),		intent(in)			::	en_k(:), ef_lst(:), T_kelvin
		real(dp),		allocatable			::	N_el_lst(:)
		integer								::	n
		!
		allocate(	N_el_lst( size(ef_lst,1))		)
		N_el_lst	=	0.0_dp
		!
		do n = 1, size(en_k)
			N_el_lst(:)	=	N_el_lst(:)		+	fd_stat(en_k(n),	ef_lst(:), T_kelvin)
		end do
		!
		return
	end function











	real(dp) pure function d0_fd_stat_deriv(e_band, e_fermi, T_kelvin)
		real(dp), 		intent(in)			::	e_band, e_fermi, T_kelvin
		real(dp)							::	T_smear, x
		!
		!	w90git:
		!if (abs (x) .le.36.0) then
        ! 	utility_w0gauss = 1.00_dp / (2.00_dp + exp ( - x) + exp ( + x) )
        !  	! in order to avoid problems for large values of x in the e
       	!else
        !	utility_w0gauss = 0.0_dp
       	!endif
       	!
       	d0_fd_stat_deriv		=	0.0_dp	
       	!
       	if(	T_kelvin > min_temp	)								then
       		T_smear			=	kBoltz_Eh_K		*	T_kelvin
       		x				=	(e_band	- e_fermi) / T_smear
       		!
       		!
       		d0_fd_stat_deriv	=	1.00_dp 	/ 	(	2.00_dp + exp( x ) + exp( -x ) 	)
       	end if

		return
	end function	

	function d1_fd_stat_deriv(fd_distrib, T_kelvin) result(fd_deriv)
		real(dp), 		intent(in)			::	fd_distrib(:,:), T_kelvin
		real(dp),		allocatable			::	fd_deriv(:,:)
		integer								::	ef_idx, n_ef, n_wf, n
		!
		!	w90git:
		!if (abs (x) .le.36.0) then
        ! 	utility_w0gauss = 1.00_dp / (2.00_dp + exp ( - x) + exp ( + x) )
        !  	! in order to avoid problems for large values of x in the e
       	!else
        !	utility_w0gauss = 0.0_dp
       	!endif
       	!		uses
       	!			2 cosh(x)	=	exp(x)	+ exp(-x)	
       	!
       	!
       	n_ef	=	size(fd_distrib,1)
       	n_wf	=	size(fd_distrib,2)
       	!
       	allocate(	fd_deriv(	n_ef,	n_wf 	))
       	fd_deriv		=	0.0_dp	
       	!
       	!
      	do n = 1, n_wf
       		do ef_idx = 1, n_ef
       			fd_deriv(ef_idx,n)	=1.00_dp 	/ 	 (		2.00_dp	+ 2.00_dp *	cosh( fd_distrib(ef_idx,n) )	)			
       		end do        		
       	end do
       	!
       	!
       	if(T_kelvin < min_temp .and. .not. t_zero_warn )	then
       		write(*,*)	"[d1_fd_stat_deriv]: WARNING fd_stat_deriv only defined for finite Temeperature (T/=0)"
       		t_zero_warn	=	.True.
       	end if
       	!
		return
	end function	











end module statistics