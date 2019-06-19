module w90_interface
#ifdef __INTEL_COMPILER
	use ifport !needed for time 
#endif
	use m_npy
	use endian_swap
	use matrix_math,					only:		is_equal_vect, is_herm_mat	
	use constants,						only:		dp, fp_acc, aUtoAngstrm, aUtoEv
	use mpi_community,					only:		mpi_id, mpi_root_id, mpi_nProcs
	use input_paras,					only:		w90_dir, 							&
													a_latt, 							&
													use_R_float, debug_mode
	implicit none
	!
	private
	public									::		read_tb_basis
	!
	character(len=64)						::		format='(a,i7.7)'
	
contains




	subroutine read_tb_basis(seed_name, R_vect, H_mat, r_mat)		!ead_tb_basis(seed_name, H_real, r_exist, r_mat)
		!	reads the following files:
		!		'_tb.dat' file			(optional)
		!		if not available
		!			'_hr.dat' file 		(mandatory)
		!			' _r.dat' file		(optional) 
		!
		character(len=*),	 			intent(in) 		::	seed_name
		real(dp),		allocatable,	intent(inout)	::	R_vect(:,:)
		complex(dp),	allocatable, 	intent(inout)	::	H_mat(:,:,:), r_mat(:,:,:,:)
		logical											::	tb_exist, hr_exist, r_exist
		!
		r_exist = .false.
		tb_exist= .false.
		hr_exist= .false.
		!
		!	FULL BASIS
		inquire(file=w90_dir//seed_name//'_tb.dat',	exist=tb_exist)
		if(tb_exist) then
			write(*,format)	'[#',mpi_id,';read_tb_basis]: "found data file '//w90_dir//seed_name//'_hr.dat'
			call read_tb_file(w90_dir//seed_name, R_vect, H_mat, r_mat)
			r_exist	= .true.
		else
			!
			!
			!	HAMILTONIAN
			inquire(file=w90_dir//seed_name//'_hr.dat', exist=hr_exist)
			if( .not. hr_exist)		then
				write(*,format)	'[#',mpi_id,';read_tb_basis]: "'//w90_dir//seed_name//'_hr.dat'//'" does not exist'
				stop '[read_tb_basis]:	ERROR could not read the real space Hamiltonian'
			end if
			call read_hr_file(w90_dir//seed_name, R_vect, H_mat)
			!
			!
			!	POSITION
			inquire(file=w90_dir//seed_name//'_r.dat', exist=r_exist)
			if(  r_exist )	then
				write(*,format)	'[#',mpi_id,'read_tb_basis]: "found data file '//w90_dir//seed_name//'_r.dat'
				call read_r_file(w90_dir//seed_name//'_r.dat', R_vect, r_mat)
			end if
			!
			!
		end if
		!
		!	CHECK IF SYMMETRYIC
		if(	.not. (size(H_mat,1) == size(H_mat,2)	)	)	stop '[read_tb_basis]:	H_mat is not symmetric'
		!
		!	CHECK BASIS CONSISTENCY
		if( r_exist	)	then
			if(size(r_mat,1)/=	3			)	stop '[read_tb_basis]:	r_mat does not live in 3D'				
			if(size(r_mat,2) /= size(H_mat,1))	stop '[read_tb_basis]:	r_mat and H_mat have different basis in first dim'
			if(size(r_mat,3) /= size(H_mat,2))	stop '[read_tb_basis]:	r_mat and H_mat have different basis in second dim'
			if(size(r_mat,4) /= size(H_mat,3))	stop '[read_tb_basis]:	r_mat and H_mat live on different real-space (cell) grid'
		end if
		!
		!	EXTENDED DEBUGGING	(HERMITICITY)
		if(debug_mode)	call real_space_debuger(R_vect, H_mat, r_mat)
		!		
		return
	end subroutine


	subroutine read_tb_file(seed_name, R_vect, tHopp, rHopp)
		character(len=*)									::	seed_name
		real(dp),		intent(out), allocatable			::	R_vect(:,:)
		complex(dp),	intent(out), allocatable			:: 	tHopp(:,:,:), rHopp(:,:,:,:)
		integer												::	mpi_unit, &
																f_nwfs, f_nSC, readLines, line, &
																cell, n
		integer												::	int3(3), index(2)
		real(dp)											::	unit_cell(3,3), compl1(2),compl3(6), rTest(3), real3(3)
		!
		mpi_unit	= 100 + mpi_id
		open(unit=mpi_unit, file=seed_name//'_tb.dat',form='formatted', status='old', action='read')
		!read unit cell
		read(mpi_unit,*)
		read(mpi_unit,*)	real3
		unit_cell(1,1:3)	= real3
		read(mpi_unit,*)	 real3
		unit_cell(2,1:3) 	= real3
		read(mpi_unit,*)	 real3
		unit_cell(3,1:3)	= real3
		unit_cell			= unit_cell / aUtoAngstrm
		!
		!overwrite the lattice read in from input
		a_latt	= unit_cell
		!read basis size
		read(mpi_unit,*)	f_nwfs
		read(mpi_unit,*)	f_nSC
		!allocate targets
		allocate(	R_vect(	3,						f_nSC		)	)
		allocate(	rHopp(	3,	f_nwfs, f_nwfs,		f_nSC		)	)
		allocate(	tHopp(		f_nwfs,	f_nwfs,		f_nSC		)	)
		!dummy read the degeneracy list
		readLines = ceiling( real(f_nSC,dp) / real(15,dp)		)
		do line = 1, readLines
			read(mpi_unit,*)
		end do
		!
		!read tHopp (eV)
		do cell = 1, f_nSC
			read(mpi_unit,*)
			!
			!	READ R-VECTOR AS ...
			if(	use_R_float	) then
				!
				!	... FLOAT
				read(mpi_unit,*)	real3(1:3)
			else
				!
				!	...	INTEGER
				read(mpi_unit,*)	int3(1:3)
				real3(:)	=		int3(:)
			end if
			R_vect(1:3,cell)	= matmul(unit_cell,	real3(1:3)	)
			!
			!	READ HOPPING
			do n = 1, f_nWfs**2
				read(mpi_unit,*)	index(1:2), compl1(1:2)
				!
				if (abs(compl1(1))	< 1e-7_dp) 	compl1(1)	=	0.0_dp
				if (abs(compl1(2))	< 1e-7_dp) 	compl1(2)	=	0.0_dp
				!
				tHopp(index(1),index(2), cell)	= cmplx(	compl1(1), compl1(2),	dp		)
			end do
		end do
		!
		!read rHopp (ang)
		do cell = 1, f_nSC
			read(mpi_unit,*)
			!
			!	read R_CELL FROM POS OP
			if(	use_R_float	) then
				read(mpi_unit,*)	real3(1:3)
			else
				read(mpi_unit,*)	int3(1:3)
				real3(:)	=		int3(1:3)
			end if
			rTest(:)	= 	matmul(unit_cell,real3(:)	)
			!
			!	CHECK IF SAME AS ORIGINAL R-CELL FROM HOPPING (see above)
			if(		 norm2( rTest(:) - R_vect(:,cell) )	> 1e-8_dp		) then
				write(*,format)	"[#",mpi_id,";read_tb_basis]: WARNING rHopp has different sc order then tHopp"
			end if
			!
			!
			do n = 1, f_nWfs**2
				read(mpi_unit,*)	index(1:2), compl3(1:6)
				!
				rHopp(1,	index(1), index(2), cell)	= cmplx( compl3(1), compl3(2),	dp	)
				rHopp(2,	index(1), index(2), cell)	= cmplx( compl3(3), compl3(4),	dp	)
				rHopp(3,	index(1), index(2), cell)	= cmplx( compl3(5), compl3(6),	dp	)
			end do
		end do 
		close(mpi_unit)
		!
		!convert to a.u.
		tHopp	= 	tHopp / aUtoEv
		rHopp	= 	rHopp / aUtoAngstrm
		!
		!
		write(*,'(a,i7.7,a,i4,a,i6,a,i6,a,i5)')	"[#",mpi_id,";read_tb_file]: unit=#",mpi_unit," SUCCESS (input interpretation: nWfs=",		&
												f_nwfs, ";	nrpts=",f_nSC,")."
		!
		!
		return
	end subroutine 




	subroutine read_hr_file(seed_name, R_vect, H_mat)
		character(len=*),	 			intent(in) 		::	seed_name
		real(dp),		allocatable,	intent(inout)	::	R_vect(:,:)
		complex(dp),	allocatable, 	intent(inout)	::	H_mat(:,:,:)
		integer											::	mpi_unit,&
															f_nwfs, f_nSC, to_read, int15(15), int3(3), m,n, idx, sc, wf
		real(dp)										::	real2(2), real3(3)
		integer,		allocatable						::	R_degneracy(:)
		!
		mpi_unit	= 	100 + mpi_id	+ mpi_nProcs

		open(unit=mpi_unit, file=seed_name//'_hr.dat',	form='formatted', action='read', access='sequential', status='old')
		!write(*,'(a,i3,a,i3)')				 		"[#",mpi_id,";read_hr_file]: 	will open _hr file on unit=",mpi_unit
		!read header
		read(mpi_unit,*)
		read(mpi_unit,*)	f_nwfs
		read(mpi_unit,*)	f_nSC
		!
		allocate(	R_degneracy(					f_nSC	)		)
		allocate(	R_vect(			3,				f_nSC	)		)
		allocate(	H_mat(			f_nWfs, f_nWfs,	f_nSC	)		)
		!
		!read wigner seitz degeneracy
		to_read = f_nSC
		idx		= 1
		do while(to_read>0)
			if(to_read >= 15)	then
				read(mpi_unit,*)		int15(1:15)
				R_degneracy(idx:idx+14)	=	int15(1:15)
				idx = idx+15	
				to_read = to_read-15
			else
				read(mpi_unit,*)		int15(1:to_read)
				R_degneracy(idx:(idx+to_read-1))	= int15(1:to_read)
				idx	= idx + to_read
				to_read	= 0
			end if
		end do
		if(	idx-1 /= f_nSC)	then
			if(mpi_id==mpi_root_id)	 then
				write(*,*)	'idx=',idx,' f_nSC=',f_nSC
				write(*,*)	R_degneracy(1:f_nSC) 
			end if


			stop '[read_hr_file]: ERROR issue reading Wigner Seitz degeneracy (idx/=f_nSC) in  _hr.dat file'
		end if
		!
		!read body
		do sc = 1, f_nSC
			do wf  = 1, f_nwfs**2
				!read next line
				if(use_R_float) then
					read(mpi_unit,*)	real3(1:3),	m, n,	real2(1:2)
				else
					read(mpi_unit,*)		int3(1:3), m,n, 	real2(1:2)
					real3	= real(int3,dp)
				end if
					!
				!get Wigner Seitz vector
				if( wf==1  ) then
					R_vect(1:3,sc)	= real3(1:3)
				else if( .not.	is_equal_vect(fp_acc,	R_vect(1:3,sc),	real3(1:3) )		) then
						write(*,*)				 		"[#",mpi_id,";read_hr_file]: 	WARNING unexpected new R_vect=",real3(:)
						write(*,'(a,i4,a,i4)',advance="no")		"	m=",m," n=",n
						write(*,*)								"	wf=",wf," sc=",sc," expected R_vect=",R_vect(:,sc)
				end if
				!
				!fill Hamiltonian
				if( abs(real2(1))< 1e-7_dp)	real2(1)	=	0.0_dp
				if( abs(real2(2))< 1e-7_dp)	real2(2)	=	0.0_dp
				H_mat(m,n, sc)	= cmplx(	real2(1),	real2(2)	, dp	)
			end do
		end do
		close(mpi_unit)
		!
		!
		!convert to a.u
		H_mat	= H_mat / aUtoEv
		!
		write(*,'(a,i7.7,a,i4,a,i6,a,i6,a,i5)')	"[#",mpi_id,";read_hr_file]: unit=#",mpi_unit," SUCCESS (input interpretation: nWfs=",		&
												f_nwfs, ";	nrpts=",size(R_vect,2),")."
		return
	end subroutine



	subroutine read_r_file(seed_name, R_vect, r_mat)
		character(len=*),	 			intent(in) 		::	seed_name
		real(dp),						intent(in)		::	R_vect(:,:)
		complex(dp),	allocatable,	intent(inout)	::	r_mat(:,:,:,:)
		integer											::	mpi_unit, &
															f_nwfs, sc, m, n, int3(3), it
		real(dp)										::	real6(6), real3(3)
		!
		mpi_unit	= 100 + mpi_id + 2*mpi_nProcs
		!
		open(unit=mpi_unit, file=seed_name//'_r.dat',	form='formatted', action='read', access='stream', status='old')
		read(mpi_unit,*)
		read(mpi_unit,*)	f_nwfs
		!
		allocate(	r_mat(3, f_nwfs, f_nwfs,	size(R_vect,2)	)			)
		!
		do sc = 1, size(R_vect,2)
			do it= 1, f_nWfs**2
				if(use_R_float) then
					read(mpi_unit,*)		real3(1:3), 		m, n, 	real6(1:6)
				else
					read(mpi_unit,*)		int3(1:3), 		m, n, 	real6(1:6)
					real3(:)	=	real(int3(:),dp)
				end if
				!
				!
				if( .not.	is_equal_vect(fp_acc,	R_vect(1:3,sc),	real3(1:3))	)	then
					write(*,*)	"[read_r_file]: R_vect=",R_vect(1:3,sc)
					write(*,*)	"[read_r_file]:	input_R=",real(int3(1:3),dp)
					stop '[read_r_file]: different R_vect order in _hr.dat and _r.dat file'
				end if
				!
				r_mat(1,m,n,sc)	=	cmplx(	real6(1)	, real6(2)	, dp	)
				r_mat(2,m,n,sc)	=	cmplx(	real6(3)	, real6(4)	, dp	)
				r_mat(3,m,n,sc)	=	cmplx(	real6(5)	, real6(6)	, dp	)
			end do
		end do
		close(mpi_unit)
		!
		!convert angstom to (au)
		r_mat	= r_mat	/ aUtoAngstrm
		!
		write(*,'(a,i7.7,a,i4,a)')	"[#",mpi_id,";read_r_file]: unit=#",mpi_unit," SUCCESS "
		!
		return
	end subroutine



!private helpers
	subroutine real_space_debuger(R_vect, H_mat, r_mat)
		real(dp),						intent(in)			::		R_vect(:,:)	
		complex(dp),	allocatable,	intent(in)			::		H_mat(:,:,:), r_mat(:,:,:,:)
		integer												::		sc, n_cells
		logical												::		hermitian
		real(dp)											::		max_err, R_tot(3)
		!
		R_tot(:) =	0.0_dp
		do sc = 1, size(R_vect,2)
			R_tot(:)	=	R_tot(:)	+	R_vect(:,sc)
		end do 
		write(*,*)	'[read_tb_basis/DEBUG-MODE]:	real space cells (input interpretation):'
		write(*,*)	"	cell	|	Rx 	,	Ry 	, 	Rz	"
		write(*,*)	"----------------------------------------------"
		do sc = 1, size(R_vect,2)
			write(*,'(i3,a,f4.1,a,f4.1,a,f4.1)')	sc,'	|		',	R_vect(1,sc),", ",R_vect(2,sc),', ',R_vect(3,sc)
		end do
		write(*,*)	"----------------------------------------------"
		write(*,'(a,f4.1,a,f4.1,a,f4.1)')	"	sum 	|	",R_tot(1),", ",R_tot(2),', ',R_tot(3)
		if(norm2(R_tot) > 1e-3_dp	) then
			write(*,*)	'[read_tb_basis/DEBUG-MODE]:	WARNING if real space cells dont sum zero hermitictiy issues might arise'
		end if
		!
		!write(*,*)	'[read_tb_basis/DEBUG-MODE]:----start debuging real space basis---------'
		n_cells	=	size(R_vect,2)
		!
		hermitian = size(H_mat,1) == size(H_mat,2)
		do sc = 1, n_cells
			!
			!			HAMILTONIAN
			!
			if( .not.	 is_herm_mat(	H_mat(  :,:,sc),max_err)) then
			 	write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space Hamiltonian is not herm in sc= ',sc,'; largest error: ', max_err
			 	hermitian	=	.false.
			 	!stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (Hamiltonian)"
			end if
			!
			!
			if(allocated(r_mat)) then
				!------------------------------------------------------------------------------------------------------------------------------------
				!			X	- POSITION																											!
				!																																	!
				if( .not.	is_herm_mat(	r_mat(1,:,:,sc),max_err)) then 																			!
					write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space x-position op. not herm, largest error: ', max_err,' sc =',sc			!
					hermitian	=	.false.																											!	
					!stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (x-pos op)"															!	
				end if																																!		
				!																																	!
				!			Y	- POSITION																											!
				!																																	!
				if( .not.	is_herm_mat(	r_mat(2,:,:,sc),max_err)) then 																			!
				write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space y-position op. not herm, largest error: ', max_err	,' sc =',sc							!
					hermitian	=	.false.																											!
					!stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (y-pos op)"															!				
				end if
				!
				!			Z	- POSITION
				!
				if( .not.	is_herm_mat(	r_mat(3,:,:,sc),max_err)) then 
					write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space z-position op. not herm, largest error: ', max_err,' sc =',sc	
					hermitian	=	.false.
					!stop "'[read_tb_basis/DEBUG-MODE]: STOP non Hermitian real space basis (z-pos op)"
				end if
				!------------------------------------------------------------------------------------------------------------------------------------
			end if
		end do
		!
		if(hermitian)	then
			write(*,*)	"[read_tb_basis/DEBUG-MODE]: SUCCESS real space basis is hermitian"
		else
			write(*,*) '[read_tb_basis/DEBUG-MODE]: ERROR real space basis IS NOT hermitian'
		end if
		!
		!
		!write(*,*)	'[read_tb_basis/DEBUG-MODE]:----finished debuging real space basis---------'
		!
		return
	end subroutine




end module w90_interface