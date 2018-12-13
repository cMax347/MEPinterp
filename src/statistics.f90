module statistics
	use constants, 		only:			dp, kBoltz_Eh_K

	implicit none

	public					::			fd_stat,		fd_stat_deriv, 			&
										fd_get_N_el, fd_count_el
	private

contains



	pure function get_fermi_energy(	N_el,	T_kelvin 	)
		integer,		intent(in)		::	N_el 		! number of electrons in system
		real(dp),		intent(in)		::	T_kelvin	! smearing paramter
		real(dp)						::	N_el_loc, N_el_glob
		real(dp)						::	e_fermi, e_fermi_max, e_fermi_min
		!
		N_el_loc	=	0.0_dp
		N_el_glob	=	0.0_dp
		!
		
		N_efermi	=	200
		e_fermi		=	e_fermi_min
		delta_e		=	e_fermi_max - e_fermi_min)	/	N_e_fermi		

		do while(e_fermi <= e_fermi_max) 

			if(	abs(	real(N_el,dp)	-	N_el_glob)	< 1e-3_dp	) then


			e_fermi	=	e_fermi + delta_e
		end do

		!
		! 	sum over local k-pts

		do kiz = 1, mp_grid(3)
			do kiy = 1, mp_grid(2)
				do kix = 1, mp_grid(1)
					ki	=	get_rel_kpt(kix, kiy, kiz, kpt)
					!
					!
					if( mpi_ki_selector(ki, num_kpts)	) then
						call get_wann_interp(do_gauge_trafo, H_tb, r_tb, a_latt, recip_latt, R_vect, ki, kpt(:), 	en_k, V_ka, A_ka, Om_kab )
						!
						N_el_loc	=	N_el_loc	+	fd_get_N_el(en_k, e_fermi, T_kelvin)
					end if
				end do
			end do
		end do
		!
		!	sum over all k-pts
		call MPI_REDUCE(	N_el_loc,  N_el_glob, 1, 	MPI_DOUBLE_PRECISION, MPI_SUM, mpi_root_id, MPI_COMM_WORLD,	ierr)






	pure subroutine fd_count_el(en_k, eFermi, T_kelvin, el_count, sum_loc, n_el_min, n_el_max)
		real(dp),		intent(in)					::		en_k(:), eFermi, T_kelvin
		real(dp),		intent(out)					::		el_count
		real(dp),		intent(inout)				::		sum_loc, n_el_min, n_el_max
		!
		el_count			=	fd_get_N_el(en_k, eFermi, T_kelvin)
		sum_loc				=	sum_loc		+	el_count
		!		!
		if(	el_count	<	n_el_min	)		n_el_min	= 	el_count
		if( el_count	>	n_el_max	)		n_el_max	=	el_count
		!
		return
	end subroutine



	real(dp) pure function fd_stat(e_band, e_fermi,	T_kelvin) 
		real(dp), 		intent(in)			::	e_band, e_fermi, T_kelvin
		real(dp)							::	T_smear
		!
		T_smear			=	kBoltz_Eh_K		*	T_kelvin
		!
		if(T_smear < 1e-6_dp ) T_smear = 1e-6_dp

		!
		fd_stat		 	= 	1.0_dp	/	(		1.0_dp	+	exp(	(e_band	- e_fermi)	/	(T_smear))				)
		!
		return
	end function




	real(dp) pure function fd_get_N_el(en_k, e_fermi, T_kelvin)
		real(dp),		intent(in)			::	en_k(:), e_fermi, T_kelvin
		integer								::	n
		!
		fd_get_N_el	=	0.0_dp
		!
		do n = 1, size(en_k)
			fd_get_N_el	=	fd_get_N_el	+	fd_stat(en_k(n),e_fermi, T_kelvin)
		end do
		!
		return
	end function



	real(dp) pure function fd_stat_deriv(e_band, e_fermi, T_kelvin)
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
       	!
       	T_smear			=	kBoltz_Eh_K		*	T_kelvin
       	x				=	(e_band	- e_fermi) / T_smear
       	!
       	!
       	fd_stat_deriv	=	1.00_dp 	/ 		(	2.00_dp + exp( x ) + exp( -x ) 			)

		return
	end function	

end module statistics