module input_paras
	use constants,					only:		dp, fp_acc, pi_dp,			&
												mpi_id, mpi_root_id, mpi_nProcs, ierr
	use matrix_math,				only:		crossP

	use mpi
	use m_config


	implicit none

	private
	public								::		&
												!routines		
												init_parameters, 	my_mkdir,						 		 		&
												get_rel_kpts, get_rel_kpt,	get_recip_latt,							&
												!dirs
												w90_dir, out_dir, raw_dir,											&
												!jobs
												plot_bands,	use_interp_kpt,											&
												!vars
												seed_name,	valence_bands,											&
												a_latt, recip_latt, mp_grid



	


	
	integer						::	valence_bands, mp_grid(3)
	character(len=3)			:: 	seed_name
	character(len=4)			::	out_dir ="out/"					
	character(len=9)			::	w90_dir	="w90files/"
	character(len=4)			::	raw_dir ="raw/"
	logical						::	plot_bands, use_interp_kpt
	real(dp)					::	a_latt(3,3), a0, recip_latt(3,3)




	contains




!public
	subroutine init_parameters()
		type(CFG_t) 			:: 	my_cfg
		real(dp)				::	a1(3), a2(3), a3(3)
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
			call CFG_add_get(my_cfg,	"unitCell%a1"      				,	a1(1:3)  	  	,	"a_x lattice vector"				)
			call CFG_add_get(my_cfg,	"unitCell%a2"      				,	a2(1:3)  	  	,	"a_y lattice vector"				)
			call CFG_add_get(my_cfg,	"unitCell%a3"      				,	a3(1:3)  	  	,	"a_z lattice vector"				)
			call CFG_add_get(my_cfg,	"unitCell%a0"					,	a0					,	"lattice scaling factor "			)
			!
			a_latt(1,1:3)	= a1(1:3)
			a_latt(2,1:3)	= a2(1:3)
			a_latt(3,1:3)	= a3(1:3)
			a_latt		=	a0 * a_latt
			![wannInterp]
			call CFG_add_get(my_cfg,	"wannInterp%mp_grid"			,	mp_grid(1:3)		,	"interpolation k-mesh"				)
			call CFG_add_get(my_cfg,	"wannInterp%seed_name"			,	seed_name			,	"seed name of the TB files			")
			call CFG_add_get(my_cfg,	"wannInterp%use_interp_kpt"		,	use_interp_kpt		,	"use w90 _geninterp.kpt file"		)
			![mep]
			call CFG_add_get(my_cfg,	"MEP%valence_bands"				,	valence_bands		,	"number of valence_bands"			)

			write(*,*)					"**********************init_parameters********************************************"
			write(*,*)					"parallelization with ",mpi_nProcs," MPI threads"
			write(*,'(a,i3,a)')			"[#",mpi_id,";init_parameters]: input interpretation:"
			write(*,*)					"[methods]"
			write(*,*)					"	plot_bands=",plot_bands
			write(*,*)					"[unitCell]"
			write(*,*)					"	a1=",a1(1:3)
			write(*,*)					"	a2=",a2(1:3)
			write(*,*)					"	a3=",a3(1:3)
			write(*,*)					"	a0=",a0
			write(*,*)					"[wannInterp]"
			write(*,*)					"	seed_name=",seed_name
			write(*,'(a,200(i3))')		"	mp_grid=",mp_grid(1:3)
			write(*,*)					"	use_interp_kpt=",use_interp_kpt
			write(*,*)					"[mep]"
			write(*,'(a,i4)')			"	val bands=",valence_bands
			

			!make the output folder
			call my_mkdir(out_dir)
			call my_mkdir(raw_dir)				
			
			write(*,*)					"---------------------------------------------------------------------------"
		end if

		!ROOT BCAST
		call MPI_BCAST(		plot_bands		,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_BCAST(		a_latt			,			9			,	MPI_DOUBLE_PRECISION	,		mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_BCAST(		use_interp_kpt	,			1			,		MPI_LOGICAL			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
		call MPI_BCAST(		valence_bands	,			1			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
		call MPI_BCAST(		seed_name(:)	,	len(seed_name)		,		MPI_CHARACTER		,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
		call MPI_BCAST(		mp_grid			,			3			,		MPI_INTEGER			,		mpi_root_id,	MPI_COMM_WORLD,	ierr)
		!
		!get reciprocal grid
		recip_latt	= get_recip_latt(a_latt)
		!
		return 
	end subroutine




	function get_recip_latt(a_latt) result(recip_latt)
		!	
		!	see eq.(2)		PRB, 13, 5188 (1976)
		!
		real(dp)				::	recip_latt(3,3)
		real(dp), intent(in)	::	a_latt(3,3)
		real(dp)				::	unit_vol, a1(3), a2(3), a3(3)
		!
		a1(1:3)		=	a_latt(1,1:3)
		a2(1:3)		=	a_latt(2,1:3)
		a3(1:3)		=	a_latt(3,1:3)
		!get Unit cell volume
		unit_vol	=	dot_product(	crossP( a1(1:3) , a2(1:3)	)		,	a3(1:3)	)	 
		!
		!
		recip_latt(1,1:3)	= 2.0_dp * pi_dp * crossP( a2(1:3) , a3(1:3) ) / unit_vol
		recip_latt(2,1:3)	= 2.0_dp * pi_dp * crossP( a3(1:3) , a1(1:3) ) / unit_vol
		recip_latt(3,1:3)	= 2.0_dp * pi_dp * crossP( a1(1:3) , a2(1:3) ) / unit_vol
		!
		return
	end function







	subroutine get_rel_kpts(mp_grid, kpt_latt)
		!get relative k-pt following the Monkhorst Pack scheme 
		!	see eq.(4)		PRB, 13, 5188 (1976)
		!
		integer,					intent(in)	::	mp_grid(3)
		real(dp),	allocatable, intent(inout)	::	kpt_latt(:,:)
		real(dp)								::	kpt(3)
		integer									::	qix, qiy, qiz, qi_idx, qi_test
		!
		allocate(	kpt_latt(3,	mp_grid(1)*mp_grid(2)*mp_grid(3)	)	)
		!
		qi_test = 0
		do qiz = 1, mp_grid(3)
			do qiy = 1, mp_grid(2)
				do qix = 1, mp_grid(1)
					qi_idx	= get_rel_kpt(qix,qiy,qiz, mp_grid, kpt	)	
					kpt_latt(:,qi_idx)	=	kpt(:)
					qi_test	= qi_test +1
					if(qi_idx /= qi_test) then
						write(*,'(a,i10,a,i10)')	'[get_rel_kpts]: WARNING k-mesh order! qi_idx=',qi_idx,' vs ',qi_test,'=qi_test'
					end if
				end do
			end do
		end do
		!
	end subroutine


	integer function get_rel_kpt(qix, qiy, qiz, mp_grid, kpt)
		integer,	intent(in)	::	qix, qiy, qiz, mp_grid(3)
		real(dp),	intent(out)	::	kpt(3)
		!
		kpt(1)	=	(	 2.0_dp*real(qix,dp)	- real(mp_grid(1),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(1),dp) )
		kpt(2)	=	(	 2.0_dp*real(qiy,dp)	- real(mp_grid(2),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(2),dp) )
		kpt(3)	=	(	 2.0_dp*real(qiz,dp)	- real(mp_grid(3),dp) - 1.0_dp		) 	/ 	( 2.0_dp*real(mp_grid(3),dp) )
		!
		!	start at 0 convention
		get_rel_kpt	=	(qix-1) + mp_grid(2)	* ( (qiy-1) + mp_grid(3) * (qiz-1)	)
		!
		!	start at 1 convention
		get_rel_kpt	=	get_rel_kpt + 1
		return
	end function


	subroutine my_mkdir(dir)
		character(len=*)			::	dir
		logical						::	dir_exists
		character(len=8)		::	mkdir="mkdir ./"	!to use with system(mkdir//$dir_path) 	
		!
		inquire(file=dir, exist=dir_exists)
		if( .not. dir_exists )	then
			call system(mkdir//dir)
			write(*,'(a,i3,a,a)')	"[#",mpi_id,"; init_parameters]: created directory ",dir
		end if
		!
		return
	end subroutine


!private




end module input_paras