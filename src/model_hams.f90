module model_hams
	use constants,		only:		dp, i_dp,	&
									aUtoEv, aUtoAngstrm
	use mpi_community,	only:		mpi_root_id, mpi_id, mpi_nProcs, ierr
	use m_config

	!
	implicit none
	!
	private
	public					::		model_ham_allocator,	get_kspace_ham
	!
	save
	integer,	parameter	::		n_free_e_gas=2,	&	
									n_rashba=2,		&
									id_rashba=0,	&
									x=1, y=2, z=3		
	!		RASHBA PARAMETER
	real(dp)					::	rashba_a, alpha_R, rashba_Vex
	real(dp)					::	rashba_nMag(3)
	!
	contains
!^
!^
!PUBLIC
	subroutine	model_ham_allocator(model_id, allo_velo,	e_k, V_ka)
		integer,						intent(in)		::	model_id
		logical,						intent(in)		::	allo_velo
		real(dp),		allocatable, 	intent(inout)	::	e_k(:)
		complex(dp),	allocatable,	intent(inout)	::	V_ka(:,:,:)
		!
		write(*,'(a,i7.7,a,i2.2)')	"[",mpi_id,";model_ham_allocator]: selected HAM #",model_id
		!
		select case(model_id)
			case(id_rashba)
				allocate(		e_k(n_rashba)		)
				call read_rashba_config(rashba_a, rashba_Vex, rashba_nMag)
			! ***PUT ADDITIONAL MODELS HERE***
			!case(id_new_case)
			!	allocate(	e_k(n_new_case)		)
			!~
			case default
				allocate(		e_k(n_free_e_gas)	)
		end select
		!
		if(	allo_velo )		&
			allocate(	V_ka(	3,	size(e_k,1), size(e_k,1)	))
		!~~~~
		return
	end subroutine
	!~
	!
	!
	subroutine get_kspace_ham(model_id,	kpt_idx, kpt,	H_k,	H_ka)
		integer,		intent(in)					::	model_id,	kpt_idx
		real(dp),		intent(in)					::	kpt(3)
		complex(dp),	allocatable,	intent(out)	::	H_k(:,:),	H_ka(:,:,:)			
		!
		select case(model_id)
			case(id_rashba)
				call	rashba_ex_model(	kpt_idx, kpt, 		H_k, H_ka)	
			! ***PUT ADDITIONAL MODELS HERE***
			!case(id_new_case)
			!	call	new_model(	kpt_idx, kpt,		H_k, H_ka)
			!~
			case default
				call	free_e_gas(			kpt_idx, kpt,		H_k, H_ka)
		end select
		!
		return
	end subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~
!~
!
!PRIVATE	-	MODEL SYSTEMS
!
	subroutine free_e_gas(	kpt_idx, kpt, H_k, H_ka)
		integer,		intent(in)						::	kpt_idx
		real(dp),		intent(in)						::	kpt(3)
		complex(dp),	allocatable, intent(inout)		::	H_k(:,:), H_ka(:,:,:)
		integer											::	ni
		complex(dp)										::	H_diago
		!
		allocate(	H_k(	n_free_e_gas,	n_free_e_gas	)	)
		allocate(	H_ka(3,	n_free_e_gas,	n_free_e_gas	)	)
		!
		!		HAMILTONIAN
		H_k			=	0.0_dp
		H_ka		=	0.0_dp
		H_diago		=	dot_product(	kpt,	kpt	)	/	2.0_dp		
		!
		do ni = 1, n_free_e_gas
			H_k(	ni,ni)	=	H_diago
			H_ka(:,	ni,ni)	=	kpt(:)	
		end do
		!	
		!
		return
	end subroutine
	!~
	!~
	!~
	subroutine read_rashba_config(aR, Vex, nMag)
		real(dp),		intent(out)	::	aR, Vex 
		real(dp),		intent(out)	::	nMag(3)
		logical						::	input_exist
		type(CFG_t) 				:: 	rash_cfg
		!
		!	INIT
		aR		=	0.0_dp
		Vex 	=	0.0_dp
		nMag	=	(/ 0.0_dp, 1.0_dp, 0.0_dp /)
		!
		!	READ CONFIG
		inquire(file="./rashba.cfg",exist=input_exist)
		if( input_exist)	then		
			!OPEN FILE
			call CFG_read_file(rash_cfg,"./rashba.cfg")
			!
			call CFG_add_get(rash_cfg,	"rashba_model%aR"	,	  aR		,	"rashba coupling strength in eV Ang"	)
			call CFG_add_get(rash_cfg,	"rashba_model%Vex"	,	 Vex		,	"exchange coupling in eV"				)
			call CFG_add_get(rash_cfg,	"rashba_model%nMag"	,  nMag(1:3)	,	"magnetization direction"				)
			!
			if(		abs(norm2(nMag))	> 1e-3_dp)	&
				nMag(:)	=	nMag(:)	/	norm2(nMag)	
			!
			write(*,'(a,i5.5,a,f6.3,a,f6.3,a,a,f3.1,a,f3.1,a,f3.1,a)')	"[#",mpi_id,";read_rashba_config]:	SUCCESS (	aR=",			&
																						aR,	" (eV Ang);	Vex=",Vex," (eV); ",		&
																				"	nmag= (/,",nMag(1)," ",nMag(2)," " ,nMag(3)," /)	)"
		else
			write(*,'(a,i5.5,a)')	"[#",mpi_id,";read_rashba_config]:	WARNING could not find ./rashba.cfg , will use standard rashba model"
		end if
		!
		!	UNIT CONVERSION
		aR		=	aR		/	(aUtoEv*aUtoAngstrm)
		Vex		=	Vex		/	aUtoEv
		!
		
		!~
		!~
		return
	end subroutine

	!
	subroutine rashba_ex_model(	kpt_idx, kpt, H_k, H_ka)
		integer,		intent(in)			::	kpt_idx
		real(dp),		intent(in)			::	kpt(3)
		complex(dp),	allocatable, intent(inout)		::	H_k(:,:), H_ka(:,:,:)
		complex(dp)										::	H_diago
		real(dp)										::	V_exxx
		integer											::	spin, xi, up=1, dw=2
		logical,		save							::	printed
		!
		!
		allocate(	H_k(	n_rashba	 ,	n_rashba	))
		allocate(	H_ka(3,	n_rashba	 ,	n_rashba	))
		!
		H_k		=	0.0_dp
		H_ka	=	0.0_dp
		H_diago	=	dot_product(	kpt(x:y),	kpt(x:y)	)	/	2.0_dp	
		!
		!^^^^^^^^^
		!
		!	KINETIC E
		do spin	= 1, n_rashba
			H_k(		spin,spin)	=	H_diago
			H_ka(x:y,	spin,spin)	=	kpt(x:y)
		end do
		!~~~~~
		!
		!	RASHBA
		H_k(	up,dw)	=	H_k(	up,dw)	+	rashba_a * (	kpt(y)	+	i_dp * kpt(x)	)	
		H_ka(x,	up,dw)	=	H_ka(x,	up,dw)	+	rashba_a *			 i_dp
		H_ka(y, up,dw)	=	H_ka(y, up,dw)	+	rashba_a 	
		!~~~~~
		!
		!	EXCHANGE
		V_exxx			=	rashba_Vex / 2.0_dp
		!
		H_k(	up,up)	=	H_k(	up,up)	+	V_exxx * 		rashba_nMag(z)
		H_k(	up,dw)	=	H_k(	up,dw)	+	V_exxx * ( rashba_nMag(x)	- i_dp * rashba_nMag(y) )
		H_k(	dw,dw)	=	H_k(	dw,dw)	-	V_exxx * 		rashba_nMag(z)
		!~~~~~
		!
		!	LOWER TRIANGLE (CC)
		H_k(dw,up)			=	conjg(		H_k(up,dw)	)
		do xi = x, y
			H_ka(xi,dw,up)	=	conjg(	H_ka(xi,up,dw)	)
		end do
		!
		return 
	end subroutine














	!~
	!~
	!~
	! ***PUT ADDITIONAL MODELS HERE***
	!
	!subroutine	new_model(	kpt_idx, kpt, H_k, H_ka)
	!	integer,		intent(in)			::	kpt_idx
	!	real(dp),		intent(in)			::	kpt(3)
	!	complex(dp),	allocatable, intent(inout)		::	H_k(:,:), H_ka(:,:)
	!			~~~~~~~
	!end subroutine
	!~
	!~
	!~
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end module 