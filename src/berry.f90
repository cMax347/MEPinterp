module berry
	!this module uses a semiclassic approach to calculate the first ordrer correction
	!	to the polariztion induced by a perturbive magnetic field
	! 	see Niu PRL 112, 166601 (2014)
	!use omp_lib
	use mpi
	use matrix_math,	only:	my_Levi_Civita
	use constants,		only:	dp, aUtoAngstrm, auToTesla,						&			
								mpi_root_id, mpi_id, mpi_nProcs, ierr			
	use statistics,		only:	fermi_dirac
	use input_paras,	only:	a_latt, 										&
								valence_bands, 									&
								seed_name,										&
								kubo_tol,										&
								hw, eFermi, T_kelvin, i_eta_smr, 				&
								unit_vol
	use file_io,		only:	mpi_read_tb_basis,								&
								write_mep_tensors,								&
								write_ahc_tensor,								&
								write_opt_tensors
	use k_space,		only:	get_recip_latt, get_mp_grid, 					&
								get_rel_kpt,									&
								normalize_k_int
	use wann_interp,	only:	get_wann_interp
	!
	use mep_niu,		only:	mep_niu_CS,	mep_niu_IC, mep_niu_LC
	use kubo,			only:	kubo_ahc_tens,									&
								kubo_opt_tens
		
	implicit none



	private
	public ::			berry_worker
	!
	save
	integer									::		num_kpts
contains



!TODO CHECK INDEXING OF VELO ACONN FCURV AFTER THEY HAVE BEEN CALCULATED


!public
	subroutine	berry_worker()
		!	interpolate the linear magnetoelectric response tensor alpha
		!
		!	lc	: local current contribution
		!	ic	: itinerant current contribution
		!	cs	: chern - simons term	
		!
		real(dp)							::	kpt(3),	recip_latt(3,3),				&
												!local sum targets:
												mep_tens_ic_loc(	3,3),	 			&
												mep_tens_lc_loc(	3,3),				&
												mep_tens_cs_loc(	3,3),				&
												kubo_ahc_loc(		3,3),				&
												!mpi sum targets:
												mep_tens_ic_glob(	3,3),				&
												mep_tens_lc_glob(	3,3),				&
												mep_tens_cs_glob(	3,3),				&
												kubo_ahc_glob(		3,3)						
												!
		complex(dp)							::	kubo_opt_s_loc(		3,3),				&
												kubo_opt_a_loc(		3,3),				&												
												kubo_opt_s_glob(	3,3),				&
												kubo_opt_a_glob(	3,3),				&
												tempS(3,3), tempA(3,3)
												!
		integer								::	kix, kiy, kiz, ki, 						&
												mp_grid(3), n_ki_loc, n_ki_glob
		complex(dp),	allocatable			::	H_tb(:,:,:), r_tb(:,:,:,:), 			&
												A_ka(:,:,:), Om_kab(:,:,:,:),			&
												V_ka(:,:,:)
		real(dp),		allocatable			::	en_k(:), R_vect(:,:)
		!
		!
		!	k-space
		mp_grid				=	get_mp_grid()
		n_ki_loc			= 	0
		num_kpts			= 	mp_grid(1)*mp_grid(2)*mp_grid(3)
		recip_latt			=	get_recip_latt()
		!
		!read real & allocate k-space
		call mpi_read_tb_basis(	seed_name, R_vect,		 	H_tb, r_tb					)
		call kspace_allocator(		H_tb, r_tb, 			en_k, V_ka, A_ka, Om_kab	)
		!
		!	MagnetoElectric Polarizability
		mep_tens_ic_loc		=	0.0_dp
		mep_tens_lc_loc		=	0.0_dp
		mep_tens_cs_loc		=	0.0_dp
		mep_tens_ic_glob	=	0.0_dp
		mep_tens_lc_glob	=	0.0_dp
		mep_tens_cs_glob	=	0.0_dp
		! 	Anomalous Hall Conducitivity
		kubo_ahc_loc		=	0.0_dp
		kubo_ahc_glob		=	0.0_dp
		!	Optical Conductivity
		kubo_opt_a_loc		=	0.0_dp
		kubo_opt_s_loc		=	0.0_dp
		kubo_opt_a_glob		=	0.0_dp
		kubo_opt_s_glob		=	0.0_dp
		!
		!	
		write(*,'(a,i3,a,i4,a)')		"[#",mpi_id,"; berry_worker]: I start interpolating now (nValence=",valence_bands,")."
		do kiz = 1, mp_grid(3)
			do kiy = 1, mp_grid(2)
				do kix = 1, mp_grid(1)
					!
					!
					ki	=	get_rel_kpt(kix, kiy, kiz, kpt)
					if( mpi_ki_selector(ki, num_kpts)	) then
						!
						!
						call get_wann_interp(H_tb, r_tb, a_latt, recip_latt, R_vect, kpt(:), 	en_k, V_ka, A_ka, Om_kab )
						!
						!MEP
						mep_tens_ic_loc	=	mep_tens_ic_loc + 	mep_niu_IC(V_ka, en_k)		!	itinerant		(Kubo)
						mep_tens_lc_loc	=	mep_tens_lc_loc + 	mep_niu_LC(V_ka, en_k)		!	local			(Kubo)
						mep_tens_cs_loc =	mep_tens_cs_loc + 	mep_niu_CS(A_ka, Om_kab)		!	chern simons	(geometrical)
						!
						!AHC
						kubo_ahc_loc	=	kubo_ahc_loc	+ 	kubo_ahc_tens(en_k,	Om_kab, eFermi, T_kelvin, unit_vol)
						!
						!OPT
						!hw, vol, e_fermi, T_kelvin, i_eta_smr , en_k, A_ka, opt_symm, opt_asym
						call kubo_opt_tens(hw, unit_vol, eFermi, T_kelvin, i_eta_smr, en_k, A_ka, 	tempS, tempA)
						kubo_opt_s_loc	=	kubo_opt_s_loc	+	tempS							
						kubo_opt_a_loc	=	kubo_opt_a_loc	+	tempA
						!
						n_ki_loc = n_ki_loc + 1
					end if
					!
					!
				end do
			end do
		end do
		write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; berry_worker]: finished interpolating ",n_ki_loc," kpts"
		!
		!
		!sum mep over global kpts
		call MPI_REDUCE(	n_ki_loc,				n_ki_glob,			1,		MPI_INTEGER,		MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	mep_tens_ic_loc,	mep_tens_ic_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		, 	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	mep_tens_lc_loc,	mep_tens_lc_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		, 	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	mep_tens_cs_loc,	mep_tens_cs_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	kubo_ahc_loc,		kubo_ahc_glob,			9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD, ierr)
		call MPI_REDUCE(	kubo_opt_s_loc,		kubo_opt_s_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD,	ierr)
		call MPI_REDUCE(	kubo_opt_a_loc,		kubo_opt_a_glob,		9,	MPI_DOUBLE_PRECISION,	MPI_SUM		,	mpi_root_id,	MPI_COMM_WORLD,	ierr)

		!
		!write some files
		if(mpi_id == mpi_root_id) 	then
			write(*,*)	"--------------------------------------------------------------------------------------------------------"	
			if(	n_ki_glob	/=	num_kpts	) stop "[berry_worker]:	ERROR n_ki_glob is not equal to the given mp_grid"
			!
			call normalize_k_int(mep_tens_ic_glob)
			call normalize_k_int(mep_tens_lc_glob)
			call normalize_k_int(mep_tens_cs_glob)
			!
			kubo_ahc_glob	=	kubo_ahc_glob		/	real(n_ki_glob,dp)
			kubo_opt_s_glob	=	kubo_opt_s_glob		/	real(n_ki_glob,dp)			
			kubo_opt_a_glob	=	kubo_opt_a_glob		/	real(n_ki_glob,dp)

			!
			call write_mep_tensors(mep_tens_ic_glob, mep_tens_lc_glob, mep_tens_cs_glob)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; berry_worker]: calculated MEP tensor on ",n_ki_glob," kpts"
			call write_ahc_tensor(kubo_ahc_glob)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; berry_worker]: calculated AHC tensor on ",n_ki_glob," kpts"
			call write_opt_tensors(kubo_opt_s_glob, kubo_opt_a_glob)
			write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; berry_worker]: calculated OPT tensor on ",n_ki_glob," kpts"
		end if
		!
		!
		return
	end subroutine









!private:
	!
	!
!HELPERS
	subroutine kspace_allocator(H_tb, r_tb, en_k, V_ka, A_ka, Om_kab)
		complex(dp),	allocatable, intent(inout)	::		H_tb(:,:,:), r_tb(:,:,:,:),					& 
															V_ka(:,:,:), A_ka(:,:,:), Om_kab(:,:,:,:)
		real(dp),		allocatable, intent(inout)	::		en_k(:)
		!
		!allocate k-space
		allocate(	en_k(						size(H_tb,2)	)	)
		allocate(	V_ka(	3,	size(H_tb,1),	size(H_tb,2)	)	)
		!
		if(	allocated(r_tb)	)	then	
			allocate(	A_ka(	  3,		size(r_tb,2),	size(r_tb,3)	)	)
			allocate(	Om_kab(	3,	3,		size(r_tb,2),	size(r_tb,3)	)	)
			write(*,'(a,i3,a)')		"[#",mpi_id,"; berry_worker]: will use position operator"	
		end if
		!
		return
	end subroutine


	logical pure function mpi_ki_selector(ki_request, num_kpts)
		integer,		intent(in)		::		ki_request, num_kpts 
		integer							::		ki_todo
		!
		mpi_ki_selector = .false.
		!
		loop_todos: do ki_todo = mpi_id +1, num_kpts, mpi_nProcs
			if(	ki_request == ki_todo) then
				mpi_ki_selector = .true.
				exit loop_todos
			end if
		end do loop_todos
		!
		return
	end function






end module berry


!old:
!	subroutine	calcFirstOrdP(polQuantum, centiMet, Bext, prefactF3, Fcurv, Aconn, Velo, En, centers_F2, centers_F3, centers_F3_essin)
!		!calculates the first order polarization p1 according to
!		!	P'= -int_dk [0.5 (Curv.Velo)*B_ext + a']
!		!
!		!	returns centers in Angstroem
!		!
!		!
!		real(dp),		intent(in)		::	polQuantum, centiMet, Bext(3), prefactF3, Fcurv(:,:,:,:), Aconn(:,:,:,:), En(:,:)		
!		complex(dp),	intent(in)		::	Velo(:,:,:,:)			
!		real(dp),		intent(out)		::  centers_F2(:,:), centers_F3(:,:), centers_F3_essin(:,:)
!		!real(dp)						::	pnF2(3), pnF3(3)
!		real(dp)						:: 	F2(3,3), F3(3,3), F2k(3,3), F3k(3,3), F3essin(3,3), F3essinK(3,3), sumF2(3), sumF3(3), &
!											p2Test(3), p3Test(3),p3_essin_test(3), p2max, p3max, p2min, p3min, kpt(3)
!		real(dp)						:: 	densCorr(3)
!		integer							:: 	n, ki, kSize, ind, k2max, k3max, k2min, k3min
!		character(len=12)				::	fname
!		!
!		kSize		= size(Velo,4)
!		!
!		if(	kSize == size(qpts,2)	)	fname = 'response.txt' 
!		if(	kSize == size(kpts,2)	)	fname = 'interpol.txt'
!		!
!		write(*,*)"[calcFirstOrdP]: start calculating P' via semiclassic approach"
!		write(*,'(a,f8.3,a,f8.3,a,f8.3,a)')"[calcFirstOrdP]: Bext=(",Bext(1)*auToTesla,", ",Bext(2)*auToTesla,", ",Bext(3)*auToTesla,") T"
!		write(*,*)"[calcFirstOrdP]: will use ",size(Velo,3)," states"
!
!		centers_F2 			= 0.0_dp
!		centers_F3 			= 0.0_dp
!		centers_F3_essin 	= 0.0_dp
!		
!		!!!!$OMP PARALLEL DEFAULT(SHARED)  &
!		!!!!$OMP PRIVATE(n, ki, densCorr, F2, F2k, F3, F3k, p2max, p2min, p3max, p3min, p2Test, p3Test)
!		!!!!$OMP DO SCHEDULE(STATIC)
!
!
!		open(unit=200,file=info_dir//'f2'//fname,action='write',status='replace')
!		open(unit=300,file=info_dir//'f3'//fname,action='write',status='replace')
!		open(unit=400,file=info_dir//'f3essin_'//fname,action='write',status='replace')
!
!		write(200,*)	"# f2 positional shift for each wf n, at each kpt"
!		write(300,*)	"# f3 positional shift for each wf n, at each kpt"
!		write(400,*)	"# f3 with the essin formalsism, flipped sign ?!"
!
!		write(200,*)	"# nStat | kpt(1:3) |	a_f2 (1:3,kpt)	(ang)"		
!		write(300,*)	"# nStat | kpt(1:3) |	a_f3 (1:3,kpt)	(ang)"
!		write(400,*)	"# nStat | kpt(1:3) |	a_f3 (1:3,kpt)	(ang)"
!		
!		write(200,*)	size(centers_F2,2)," ",kSize 
!		write(200,*)	Bext(3)*auToTesla
!		write(300,*)	size(centers_F3,2)," ",kSize
!		write(300,*)	Bext(3)*auToTesla
!		write(400,*)	size(centers_F3_essin,2)," ",kSize
!		write(400,*)	Bext(3)*auToTesla
!
!
!		do n = 1, size(centers_F2,2)
!			F2 = 0.0_dp
!			F3 = 0.0_dp
!			F3essin	= 0.0_dp
!			!
!			k2max = -1
!			k3max = -1
!			k2min = -1
!			k3min = -1
!
!			!GET RESPONSE MATRIX
!			p2max	= 0.0_dp
!			p3max	= 0.0_dp
!			p2min	= 100.0_dp
!			p3min	= 100.0_dp
!
!			do ki = 1, kSize		
!				!
!				!PHASE SPACE DENSITY CORRECTION
!				densCorr(1:3)	= 0.5_dp * dot_product(		Fcurv(1:3,n,n,ki), Aconn(1:3,n,n,ki)	)		* Bext
!				if( norm2(densCorr) > acc ) write(*,*)	"[calcFirstOrdP]: WARNING the densCorr is none zero, norm2(densCorr)=",norm2(densCorr)
!				!
!				!POSITIONAL SHIFT
!				call getF2(n,ki,Velo,En, F2k)
!				call getF3(prefactF3, n,ki,Velo,En, F3k)
!				call getF3essin(prefactF3, n, ki, Velo, En, F3essinK)
!				!sum over K
!				F2 = F2 + F2k
!				F3 = F3 + F3k
!				F3essin = F3essin + F3essinK
!				!
!				!search for extrema
!				p2Test = matmul(F2k,Bext)* aUtoAngstrm
!				p3Test = matmul(F3k,Bext)* aUtoAngstrm
!				p3_essin_test = matmul(F3essinK,Bext) * aUtoAngstrm
!
!				kpt = 0.0_dp
!				if( kSize == size(qpts,2)	)	kpt(1:2) 	= qpts(1:2,ki)
!				if( kSize == size(kpts,2)	)	kpt(1:2)	= kpts(1:2,ki)
!				
!				write(200,'(i3,a,i5,a,e16.9,a,e16.9,a,e16.9,a,e16.9,a,e16.9,a,e16.9)')	n," ",ki," ",kpt(1)," ",kpt(2)," ",kpt(3)," ",p2Test(1)," ",p2Test(2)," ",p2Test(3)
!				write(300,'(i3,a,i5,a,e16.9,a,e16.9,a,e16.9,a,e16.9,a,e16.9,a,e16.9)')	n," ",ki," ",kpt(1)," ",kpt(2)," ",kpt(3)," ",p3Test(1)," ",p3Test(2)," ",p3Test(3)
!				write(400,'(i3,a,i5,a,e16.9,a,e16.9,a,e16.9,a,e16.9,a,e16.9,a,e16.9)')	n," ",ki," ",kpt(1)," ",kpt(2)," ",kpt(3)," ",p3Test(1)," ",p3Test(2)," ",p3_essin_test(3)
!
!				if( norm2(p2Test) > p2max) then
!					p2max = norm2(p2Test)
!					k2max = ki
!				else if (norm2(p2Test) < p2min) then
!					p2min = norm2(p2Test)
!					k2min = ki
!				end if
!				if( norm2(p3Test) > p3max) then
!					p3max = norm2(p3Test)
!					k3max = ki
!				else if (norm2(p3Test) < p3min) then
!					p3min = norm2(p3Test)
!					k3min = ki
!				end if
!			end do
!			!
!			write(*,'(a,i2,a,i5,a,e13.4)')	"[calcFirstOrdP]: n=",n," largest F2 shift (at #kpt=",k2max,"): ",p2max, "(ang)"
!			write(*,'(a,i2,a,i5,a,e13.4)')	"[calcFirstOrdP]: n=",n," smalles F2 shift (at #kpt=",k2min,"): ",p2min, "(ang)"
!			write(*,'(a,i2,a,i5,a,e13.4)')	"[calcFirstOrdP]: n=",n," largest F3 shift (at #kpt=",k3max,"): ",p3max, "(ang)"
!			write(*,'(a,i2,a,i5,a,e13.4)')	"[calcFirstOrdP]: n=",n," smalles F3 shift (at #kpt=",k3min,"): ",p3min, "(ang)"
!			!NORMALIZE
!			F2 = F2 / real(kSize,dp)
!			F3 = F3  / real(kSize,dp)
!			F3essin = F3essin / real(ksize,dp) 
!		
!			!APPLY MATRIX 
!			centers_F2(:,n) = matmul(F2,Bext) * aUtoAngstrm
!			centers_F3(:,n) = matmul(F3,Bext) * aUtoAngstrm
!			centers_F3_essin(:,n) = matmul(F3essin,Bext) * aUtoAngstrm
!			!
!		end do
!
!
!		close(200)
!		close(300)
!		close(400)
!
!
!		!!!!$OMP END DO
!		!!!!$OMP END PARALLEL
!		!
!
!		do ind = 1, 3
!			sumF2(ind) 	= sum( 	centers_F2(ind,:)		)  	* polQuantum *centiMet
!			sumF3(ind)	= sum(	centers_F3(ind,:) 		)	* polQuantum *centiMet		
!		end do
!		!
!		!PRINT F2
!		write(*,*)															"[calcFirstOrdP]: F2 matrix contribution:"
!		write(*,*)															" #state | 		<r>[Å]			| 		p[mu C / cm]"
!		do n = 1, size(centers_F2,2)	
!			write(*,'(i3,a,e13.4,a,e13.4,a,e13.4,a,a,e13.4,a,e13.4,a)')		n," | ", centers_F2(1,n),", ",centers_F2(2,n), ", ", centers_F2(3,n)," | ", &
!																					" (",	centers_F2(1,n) * polQuantum * centiMet		,&
!																					", ",	centers_F2(2,n) * polQuantum * centiMet		,")"
!		end do
!		write(*,'(a,e13.4,a,e13.4,a,e13.4,a)')								"sum | 						|	(", sumF2(1),", ",sumF2(2), ", ", sumF2(3),")."
!		!
!		!PRINT F3
!		write(*,*)															"[calcFirstOrdP]: F3 matrix contribution:"
!		write(*,*)															" #state | 		<r>[Å]			| 		p[mu C / cm]"
!		do n = 1, size(centers_F3,2)	
!			write(*,'(i3,a,e13.4,a,e13.4,a,e13.4,a,a,e13.4,a,e13.4,a)')		n," | ", centers_F3(1,n),", ",centers_F3(2,n),", ", centers_F3(3,n)," | ", &
!																					" (",	centers_F3(1,n)*aUtoAngstrm * polQuantum * centiMet		, &
!																					", ",	centers_F3(2,n)*aUtoAngstrm * polQuantum * centiMet		,")"
!		end do
!		write(*,'(a,e13.4,a,e13.4,a,e13.4,a)')								"sum | 						|	(", sumF3(1),", ",sumF3(2), ", ", sumF3(3),")."
!		!
!		!PRINT TOT
!		write(*,*)															"[calcFirstOrdP] total first order pol:"
!		write(*,'(a,e13.4,a,e13.4,a,e13.4,a)')								"p'= (",sumF2(1)+sumF3(1),", ",sumF2(2)+sumF3(2),", ",sumF2(3)+sumF3(3),") [mu C/ cm] "
!
!		!
!		!DEBUG
!		if(	kSize /= size(En,2)	)	 write(*,*)"[calcFirstOrdP]: WARNING Energy and velocities live on different k meshes!"
!		!
!		!
!		return
!	end subroutine
!
!	subroutine	calcFmat(prefactF3, nZero,ki, Velo ,En, Fmat)
!		!calculates the linear response F matrix for magnetic field perturbation
!		!F is derived in the semiclassical wavepacket approach (again see Niu PRL 112, 166601 (2014))
!		integer,		intent(in)		:: nZero, ki
!		complex(dp),	intent(in)		:: Velo(:,:,:,:)  !V(3,nWfs,nWfs,nK)
!		real(dp),		intent(in)		:: prefactF3, En(:,:)			!En(nK nWfs)
!		real(dp),		intent(out)		:: Fmat(3,3)
!		real(dp)						:: F2(3,3), F3(3,3)
!		!
!		Fmat	=	0.0_dp		
!		call getF2(nZero, ki, Velo, En, F2)
!		call getF3(prefactF3, nZero, ki, Velo, En, F3)
!		!
!		!
!		Fmat	=	F2 + F3
!		!
!		return
!	end subroutine
!
!
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

