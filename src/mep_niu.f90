module mep_niu
	!this module uses a semiclassic approach to calculate the first ordrer correction
	!	to the polariztion induced by a perturbive magnetic field
	! 	see Niu PRL 112, 166601 (2014)
	!use omp_lib
	use mpi
	use parameters,		only:	myLeviCivita,									&
								dp, aUtoAngstrm, auToTesla,	machineP,			&
								mpi_root_id, mpi_id, mpi_nProcs, ierr,			&
								seed_name,										&
								plot_bands,										&
								a_latt, valence_bands			
	use file_io,		only:	read_k_mesh, read_tb_basis,						&
								write_en_binary, read_en_binary,				&
								write_en_global,								&
								write_mep_tensor
	use wann_interp,	only:	wann_interp_ft, velo_interp
	use matrix_math,	only:	zheevd_wrapper

	implicit none



	private
	public ::			mep_interp

	real(dp),		parameter 	::	elemCharge	 	= 1.6021766208 * 1e-19_dp  *1e+6_dp! in  mu Coulomb

	integer									::		num_kpts
	real(dp),		dimension(3,3)			::		real_lattice, recip_lattice
	real(dp),		allocatable				::		kpt_latt(:,:)


	contains



!TODO CHECK INDEXING OF VELO ACONN FCURV AFTER THEY HAVE BEEN CALCULATED


!public
	subroutine	mep_interp()
		real(dp)							::	F2(3,3), F3(3,3),						&
												mep_tens_k(		3,3),	 				&
												mep_tens_loc(	3,3),					& 
												mep_tens_glob(	3,3) 	
		integer								::	ki, ki_loc_count
		complex(dp),	allocatable			::	H_real(:,:,:), r_mat(:,:,:,:), 			&
												U_k(:,:), H_ka(:,:,:), A_ka(:,:,:), 	&
												v_k(:,:,:)
		real(dp),		allocatable			::	en_k(:), R_vect(:,:)

		!
		!	
		!get (realtive) interp mesh
		call read_k_mesh(seed_name, kpt_latt)
		num_kpts	= 	size(kpt_latt,2)
		!
		!get real space matrice(s)
		call read_tb_basis(seed_name, R_vect, H_real, r_mat)
		!
		!allocate k-space
		allocate(	U_k(		size(H_real,1),	size(H_real,2)	)	)
		allocate(	en_k(					size(H_real,2)		)	)
		allocate(	H_ka(	3,	size(H_real,1),	size(H_real,2)	)	)
		allocate(	v_k(	3,	size(H_real,1),	size(H_real,2)	)	)
		if(	allocated(r_mat)	)		allocate(	A_ka(	3,	size(r_mat,2),	size(r_mat,3)	)		)
		!
		!
		!loop kpts
		mep_tens_loc	=	0.0_dp
		mep_tens_glob	=	0.0_dp
		ki_loc_count	=	0
		write(*,'(a,i3,a)')		"[#",mpi_id,"; mep_interp]: I start interpolating now"
		do ki = mpi_id + 1, num_kpts,	mpi_nProcs
			!
			!
			!interpolate
			call wann_interp_ft(H_real, r_mat,  R_vect,	 kpt_latt(1:3,ki),	U_k,	H_ka, A_ka)	
			!get energies
			call zheevd_wrapper(U_k, en_k)
			!get velos
			call velo_interp(U_k, en_k, H_ka, A_ka, V_k)
			!get MEP_tensor
			call get_F2(V_k, en_k, F2)
			call get_F3(V_k, en_k, F3)
			mep_tens_k	=	F2 + F3			!
			!sum MEP over kpts
			mep_tens_loc = mep_tens_loc + mep_tens_k
			!
			!
			!write some files
			if(plot_bands)	call write_en_binary(ki, en_k)
			write(*,'(a,i3,a,i8)')	"[#",mpi_id,"; mep_interp]: interpolated ki=",ki
			ki_loc_count = ki_loc_count + 1
		end do
		write(*,'(a,i3,a,i8,a)')		"[#",mpi_id,"; mep_interp]: finished interpolating ",ki_loc_count," kpts"
		!

		!sum over all kpts
		call MPI_REDUCE(	mep_tens_loc, mep_tens_glob,	 9,	MPI_DOUBLE_PRECISION,	MPI_SUM	, 	mpi_root_id,	MPI_COMM_WORLD, ierr)
		!
		!write some files
		if(mpi_id == mpi_root_id) 	then
			write(*,*)	"---------------------------------------------------------------------------------------------------------------------------"
			call write_mep_tensor(mep_tens_glob)			
			!
			if(plot_bands)	then
				call write_en_global(kpt_latt)
			else
				!only write mep tensor if a full BZ calculation is done
				call write_mep_tensor(mep_tens_glob)
			end if
		end if
		!
		return
	end subroutine












!private:
	subroutine get_F2(velo, En,	F2)
		complex(dp),		intent(in)		::	velo(:,:,:)
		real(dp),			intent(in)		::	En(:)
		real(dp),			intent(out)		::	F2(3,3)
		integer								::	n0, m, n, 		&
												i, j, k, l
		complex(dp)							::	velo_nom
		real(dp) 							::	en_denom
		!
		F2	= 0.0_dp
		do n0 = 1, valence_bands
			!
			!MIXING
			do m = 1, size(velo,2)
				do n = 1, size(velo,2)
				 	if( n/= n0 .and. m/=n0 )	then
				 		!
				 		!TRIPLE PRODUCT
				 		do j = 1, 3
				 			do i = 1, 3
				 				do l = 1, 3
				 					do k = 1, 3
				 						velo_nom	=	velo(i,n0,n) * velo(k,n,m) * velo(l,m,n0)
				 						!
				 						en_denom	=	(	en(n0) - en(n)	)**2		 * 		(	en(n0) - en(m)	)  
				 						!
				 						if( abs(en_denom) > 1e-15_dp) then
				 							F2(i,j)	= F2(i,j) - real(myLeviCivita(j,k,l),dp) * dreal(	velo_nom	 )	/	en_denom
				 						else
				 							write(*,*)	"[get_F2]: degenerate band warning"
				 						end if	
				 					end do
				 				end do
				 			end do
				 		end do
				 		!
				 		!
				 	end if
				end do
			end do
			!
			!
		end do
		return
	end subroutine




	subroutine get_F3(velo, En,	F3)
		complex(dp),		intent(in)		::	velo(:,:,:)
		real(dp),			intent(in)		::	En(:)
		real(dp),			intent(out)		::	F3(3,3)
		integer								::	n0, n, 		&
												i, j, k, l
		complex(dp)							::	velo_nom
		real(dp) 							::	en_denom
		!
		F3	= 0.0_dp
		do n0 = 1, valence_bands
			!
			!MIXING
			do n = 1, size(velo,2)
			 	if( n/= n0 )	then
					!
					!TRIPLE PRODUCT			 		
			 		do j = 1, 3
			 			do i = 1, 3
			 				do l = 1, 3
			 					do k = 1, 3
			 						velo_nom	=	velo(i,n0,n) * velo(k,n,n0) * velo(l,n0,n0)
			 						!
			 						en_denom	=	(	en(n0) - en(n)	)**3		
			 						!
			 						if(abs(en_denom) > 1e-15_dp) then 
			 							F3(i,j)	= F3(i,j) + real(myLeviCivita(j,k,l),dp) * dreal(	velo_nom	 )	/	en_denom	
			 						else
			 							write(*,*)	"[get_F3]: degenerate band warning"
			 						end if
			 					end do
			 				end do
			 			end do
			 		end do
			 	end if
			 	!
			 	!
			end do
			!
			!
		end do
		return
	end subroutine




end module mep_Niu


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
!										F2(i,j) 	= F2(i,j) +   real(myLeviCivita(j,k,l),dp) *  dreal( Vtmp )  / eDiff	
!										!F2(i,j) 	= F2(i,j) +   myLeviCivita(j,k,l) *  dreal( Vtmp )  / eDiff	
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
!								F3(i,j) 	= F3(i,j) + real(prefactF3,dp) * real(myLeviCivita(j,k,l),dp) *	 dreal( Vtmp ) / eDiff
!								!F3(i,j) 	= F3(i,j) + prefactF3 * myLeviCivita(j,k,l) *	 dreal( Vtmp ) / eDiff
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
!								F3(i,j) 	= F3(i,j) + real(prefactF3,dp) * real(myLeviCivita(j,k,l),dp) *	 dreal( Vtmp ) / eDiff
!								!F3(i,j) 	= F3(i,j) + prefactF3 * myLeviCivita(j,k,l) *	 dreal( Vtmp ) / eDiff
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

