module parameters
	use mpi
	use m_config
	use matrix_math,				only:		crossP


	implicit none

	private
	public								::		&
												!routines		
												myExp, myLeviCivita, init_parameters, get_rel_kpts,					&
												!dirs
												w90_dir, out_dir,													&
												!jobs
												plot_bands,	use_interp_kpt,											&
												!constatns
												dp, machineP, acc,													&
												PI_dp, i_dp,														&
												aUtoAngstrm, aUtoEv, aUtoTesla, speedOfLight,						&
												!mpi vars
												mpi_root_id, mpi_id, mpi_nProcs, ierr,								&
												!vars
												seed_name,	do_gauge_trafo,	 valence_bands,							&
												a_latt, unit_vol, recip_latt, mp_grid


	!for clean double precision convention through the code
	integer, 		parameter 	:: 	dp 				= kind(0.d0)
	real(dp), 		parameter	:: 	machineP 		= 1e-15_dp
	real(dp)					:: 	acc				= 1e-14_dp
	
	!mathematical constants
	real(dp), 		parameter 	::	PI_dp 			= 4 * datan (1.0_dp)
	complex(dp),	parameter 	::	i_dp 			= dcmplx(0.0_dp, 1.0_dp)


	!physical constants
	real(dp),		parameter 	::	aUtoAngstrm 	= 0.52917721092_dp
	real(dp),		parameter 	::	aUtoEv	 		= 27.211385_dp
	real(dp),		parameter	::	aUtoTesla		= 235051.76_dp
	real(dp),		parameter	::	speedOfLight	= 137.035999_dp !in atomic units


	!
	integer						::	mpi_id, mpi_root_id, mpi_nProcs, ierr,												&
									valence_bands, mp_grid(3)
	character(len=10)			:: 	seed_name
	character(len=9)			::	w90_dir	="w90files/"
	character(len=7)			::	out_dir ="MEPout/"
	logical						::	r_exists, do_gauge_trafo, plot_bands, use_interp_kpt
	real(dp)					::	a_latt(3,3), a0, unit_vol, recip_latt(3,3)




	contains




!public
	subroutine init_parameters()
		real(dp)				::	tmp(3), a1(3), a2(3), a3(3)
		type(CFG_t) :: 	my_cfg
		!
		!ROOT READ
		if(mpi_id == mpi_root_id) then
			!OPEN FILE
			call CFG_read_file(my_cfg,"./input.txt")
			!
			![methods]
			call CFG_add_get(my_cfg,	"jobs%plot_bands"				,	plot_bands			,	"if true do a bandstructure run"	)
			!READ SCALARS
			![unitCell]
			call CFG_add_get(my_cfg,	"unitCell%a1"      				,	a_latt(1,1:3)  	  	,	"a_x lattice vector"				)
			call CFG_add_get(my_cfg,	"unitCell%a2"      				,	a_latt(2,1:3)  	  	,	"a_y lattice vector"				)
			call CFG_add_get(my_cfg,	"unitCell%a3"      				,	a_latt(3,1:3)  	  	,	"a_z lattice vector"				)
			call CFG_add_get(my_cfg,	"unitCell%a0"					,	a0					,	"lattice scaling factor "			)
			!		
			a_latt		=	a0 * a_latt
			![wannInterp]
			call CFG_add_get(my_cfg,	"wannInterp%do_gauge_trafo"		,	do_gauge_trafo		,	"if true calc velos in Ham gauge"	)
			call CFG_add_get(my_cfg,	"wannInterp%use_interp_kpt"		,	use_interp_kpt		,	"use w90 _geninterp.kpt file"		)
			call CFG_add_get(my_cfg,	"wannInterp%seed_name"			,	seed_name			,	"seed name of the Tight Binding files")
			![mep]
			call CFG_add_get(my_cfg,	"MEP%valence_bands"				,	valence_bands		,	"number of valence_bands"			)
		end if

		!ROOT BCAST
		call MPI_BCAST(		plot_bands		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_BCAST(		a_latt			,			9			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_BCAST(		do_gauge_trafo	,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
		call MPI_BCAST(		valence_bands	,			1			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
		call MPI_BCAST(		seed_name		,	len(seed_name)		,		MPI_CHAR			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
	
		a1(1:3)		=	a_latt(1,1:3)
		a2(1:3)		=	a_latt(2,1:3)
		a3(1:3)		=	a_latt(3,1:3)
		!get Unit cell volume
		tmp(1:3)	= 	crossP( a1(1:3) , a2(1:3)	)
		unit_vol	=	dot_product(	tmp(1:3)		,	a3(1:3)	)	 
		!
		!get reciprocal grid
		recip_latt(1,1:3)	= 2.0_dp * pi_dp * crossP( a2(1:3) , a3(1:3) ) / unit_vol
		recip_latt(2,1:3)	= 2.0_dp * pi_dp * crossP( a3(1:3) , a1(1:3) ) / unit_vol
		recip_latt(3,1:3)	= 2.0_dp * pi_dp * crossP( a1(1:3) , a2(1:3) ) / unit_vol
		!
		return 
	end subroutine








	complex(dp) function myExp(x)
		!supposed to boost performance
		real(dp), intent(in) :: x
		!
		myExp = dcmplx(  dcos(x) , dsin(x) ) 
		!
		return
	end function




	integer function myLeviCivita(i,j,k)
		!Hard coded Levi Civita tensor
		integer,		intent(in)		:: i,j,k
		logical							:: even, odd
		!
		!
		even	= (i==1 .and. j==2 .and. k==3) .or. (i==2 .and. j==3 .and. k==1) .or. (i==3 .and. j==1 .and. k==2)
		odd		= (i==3 .and. j==2 .and. k==1) .or. (i==1 .and. j==3 .and. k==2) .or. (i==2 .and. j==1 .and. k==3)
		!
		if(even) then
			myLeviCivita	=  1
		else if (odd) then
			myLeviCivita	= -1
		else
			myLeviCivita	=  0
		end if
		!
		!DEBUGGING
		if(even .and. odd) then
			write(*,*)"[myLeviCivita]: myLeviCivita detected even and odd, contact the programer he fucked up"
		end if
		!
		return
	end function



	subroutine get_rel_kpts(mp_grid, kpt_latt)
		!get relative k-pt following the Monkhorst Pack scheme PRB, 13, 5188 (1976)
		integer,					intent(in)	::	mp_grid(3)
		real(dp),	allocatable, intent(inout)	::	kpt_latt(:,:)
		integer									::	qix, qiy, qiz, qi_idx
		!
		allocate(	kpt_latt(3,	mp_grid(1)*mp_grid(2)*mp_grid(3)	)	)
		!
		qi_idx = 1
		do qiz = 1, mp_grid(3)
			do qiy = 1, mp_grid(2)
				do qix = 1, mp_grid(1)	
					kpt_latt(1,qi_idx)	=	(	 2.0_dp*qix	- mp_grid(1) - 1.0_dp		) 	/ 	( 2.0_dp*mp_grid(1) )
					kpt_latt(2,qi_idx)	=	(	 2.0_dp*qiy	- mp_grid(2) - 1.0_dp		) 	/ 	( 2.0_dp*mp_grid(2) )
					kpt_latt(3,qi_idx)	=	(	 2.0_dp*qiz	- mp_grid(3) - 1.0_dp		) 	/ 	( 2.0_dp*mp_grid(3) )
					qi_idx				= 	qi_idx + 1
				end do
			end do
		end do
		!
	end subroutine






end module