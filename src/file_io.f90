module file_io
	use parameters,						only:		dp, aUtoAngstrm, aUtoEv, 			&
													mpi_id, mpi_root_id,				&
													w90_dir, out_dir, raw_dir, 			&
													a_latt, mp_grid, 					&
													plot_bands,  use_interp_kpt,		&
													get_rel_kpts
	use matrix_math,					only:		is_equal_vect												

	implicit none

	private
	public									::		read_k_mesh, read_tb_basis,			& 
													write_en_binary, read_en_binary,	&
													write_en_global,					&
													write_mep_tensors

	character(len=64)						::		format='(a,i7.7)'
	integer									::		num_bands
	contains






!public read:
	subroutine read_k_mesh(seed_name,	kpt_latt)
		!	reads the file seed_name//.nnkp 
		character(len=*), 					intent(in) 			:: 	seed_name
		real(dp),			allocatable,	intent(inout)		::	kpt_latt(:,:)
		logical													::	found_kpts_pl, found_kpts_wann
		!
		found_kpts_pl	= .false.
		found_kpts_wann	= .false.
		!
		!TRY TO READ EXISTING KPT FILE
		if( plot_bands	) then
			found_kpts_pl = read_kptsgen_pl_file(kpt_latt)
		else if( use_interp_kpt ) then
			found_kpts_wann = read_geninterp_kpt_file(seed_name, kpt_latt)
		end if
		!
		!
		!IF NOT FOUND SET UP MP GRID
		if( .not. plot_bands .and. .not. found_kpts_wann) then
			if( use_interp_kpt)	then
				write(*,'(a,i3,a)',advance="no")	"[#",mpi_id,"; read_k_mesh]: no kpt file found, will generate ("
				write(*,'(i3,a,i3,a,i3,a)')									mp_grid(1),"x",mp_grid(2),"x",mp_grid(3),") MP grid"
			else 
				write(*,'(a,i3,a)',advance="no")	"[#",mpi_id,"; read_k_mesh]: a new kpt mesh will be generated ("
				write(*,'(i3,a,i3,a,i3,a)')							mp_grid(1),"x",mp_grid(2),"x",mp_grid(3),") MP grid"
			end if
			call get_rel_kpts(mp_grid, kpt_latt)
			!
			!write to file
			if( mpi_id == mpi_root_id)	call write_geninterp_kpt_file(seed_name, kpt_latt)
		end if
		!
		return
	end subroutine







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
		!inquire(file=)
		inquire(file=w90_dir//seed_name//'_tb.dat',	exist=tb_exist)
		if(tb_exist) then
			call read_tb_file(w90_dir//seed_name, R_vect, H_mat, r_mat)
			r_exist	= .true.
		else
			inquire(file=w90_dir//seed_name//'_hr.dat', exist=hr_exist)
			if( .not. hr_exist)		stop 'ERROR could not read the real space Hamiltonian'
			call read_hr_file(w90_dir//seed_name, R_vect, H_mat)
			!
			inquire(file=w90_dir//seed_name//'_r.dat', exist=r_exist)
			if(  r_exist )	then
				call	read_r_file(w90_dir//seed_name, R_vect, r_mat)
			end if
		end if
		!
		num_bands=size(H_mat,1)
		!
		!
		return
	end subroutine


	subroutine read_en_binary(qi_idx, e_bands)
		integer,		intent(in)		::	qi_idx
		real(dp),		intent(out)		::	e_bands(:)
		character(len=24)				::	filename
		!
		write(filename, format) raw_dir//'/enK.',qi_idx
		open(unit=211, file = filename, form='unformatted', action='read', access='stream',		status='old'		)
		read(211)	e_bands(:)
		close(211)
		!
		return
	end subroutine
	





!public write:
	subroutine write_geninterp_kpt_file(seed_name, kpt_latt)
		character(len=*),	intent(in)			::	seed_name
		real(dp),			intent(in)			::	kpt_latt(:,:)	
		integer									::	qi_idx, x		


		open(unit=200, file = seed_name//'_geninterp.kpt', form='formatted', action='write',	access='stream', status='replace')
		write(200,*)	'# fractional kpts created by MEPInterp'
		write(200,*)	'frac'
		write(200,*)	size(kpt_latt,2)


		do qi_idx = 1, size(kpt_latt,2)
			write(200,'(i3,a)',advance="no")	qi_idx, ' '
			write(200,'(*(f16.8))')		(kpt_latt(x,qi_idx), x= 1, size(kpt_latt,1))
		end do
		close(200)
		write(*,'(a,i3,a,a)')		'[#',mpi_id,'; write_geninterp_kpt_file]: wrote ',seed_name//'_geninterp.kpt file'

		return
	end subroutine



	subroutine write_en_binary(qi_idx, e_bands)
		integer,		intent(in)		::	qi_idx
		real(dp),		intent(in)		::	e_bands(:)
		character(len=24)				::	filename
		!
		write(filename, format) raw_dir//'/enK.',qi_idx
		open(unit=210,	file = filename, form='unformatted', action='write', access='stream',	status='replace'		)
		write(210)	e_bands(:)
		close(210) 
		write(*,'(a,i3,a,a)')	"[#",mpi_id," ;write_en_binary]: prepare bands, wrote binary file ",filename
		!
		return
	end subroutine



	subroutine write_en_global(kpt_latt)
		real(dp),	intent(in)			::	kpt_latt(:,:)
		real(dp),	allocatable			::	ek_bands(:)
		integer							::	qi_idx, band, x
		!
		allocate(	ek_bands(num_bands)		)
		!
		open(unit=220, file=out_dir//'eBands.dat', form='formatted', action='write', access='stream', status='replace')
		write(220,*)	'# energies interpolated by MEPinterp program'
		write(220,*)	'# first 3 columns give the relative k-mesh, 4th column are the enegies'
		write(220,*)	'# Kpt_idx  K_x (frac)       K_y (frac)        K_z (frac)       Energy (Hartree)  '
		do qi_idx = 1, size(kpt_latt,2)
			call read_en_binary(qi_idx, ek_bands)
			!
			do band = 1, size(ek_bands,1)
				write(220,'(200(f16.8))',advance="no")	(kpt_latt(x,qi_idx), x=1,size(kpt_latt,1)	)
				write(220, '(a,f18.8)')			'	',ek_bands(band)
			end do
			!		
		end do
		close(220)
		write(*,'(a,i3,a)')		"[#",mpi_id,"; write_en_global]: success!"
		!
		return
	end subroutine




	subroutine write_mep_tensors(mep_ic, mep_lc, mep_cs)
		real(dp),			intent(in)		::	mep_ic(3,3), mep_lc(3,3), mep_cs(3,3)
		real(dp)							::	mep_tot(3,3)
		character(len=12)					::	fname
		character(len=50)					::	info_string
		!
		!-------------------------------------itinerant contribution MEP tensor-------------
		fname		= 'mep_ic.dat'
		info_string	= '# itinerant contribution of mep tensor'
		!
		call	write_mep_file(fname,	mep_ic,	info_string )
		!
		!
		!-------------------------------------local contribution MEP tensor-----------------
		fname		= 'mep_lc.dat'
		info_string	= '# local contribution of mep tensor'
		!
		call	write_mep_file(fname,	mep_lc,	info_string )
		!
		!
		!-------------------------------------Chern-Simons term MEP tensor------------------
		fname		= 'mep_cs.dat'
		info_string	= '# Chern-Simons term of mep tensor'
		!
		call	write_mep_file(fname,	mep_cs,	info_string )
		!
		!
		!-------------------------------------total MEP tensor------------------------------
		fname		= 'mep_tens.dat'
		info_string	= '# total mep tensor (mep_tot= mep_ic+mep_lc+mep_cs)'
		!
		mep_tot		= mep_ic +	mep_lc	+	mep_cs
		call	write_mep_file(fname,	mep_tot,	info_string )
		!-----------------------------------------------------------------------------------
		!
		return
	end subroutine



!private write
	subroutine write_mep_file(fname,mep_tens, info_string)
		character(len=*), 	intent(in)		::	fname, info_string
		real(dp),			intent(in)		::	mep_tens(3,3)
		integer								::	row, clm		!
		!write result to file
		open(unit=250, file=out_dir//fname, form='formatted', 	action='write', access='stream',	status='replace')
			write(250,*)	info_string
			write(250,*)	'begin mep'
			do row = 1, 3
				write(250,'(200(f16.8,a))')		(		mep_tens(row,clm), ' ', clm=1,3)
			end do
			write(250,*)	'end mep'
		close(250)
		write(*,'(a,i3,a)')	"[#",mpi_id,"; write_mep_tensor]: success!"
		!
		!
		return
	end subroutine

























!private read:	
	logical function read_kptsgen_pl_file(kpt_latt)
		!reads the kpt file generatred by kptsgen.pl
		real(dp),			allocatable,	intent(inout)		::	kpt_latt(:,:)
		character(len=50)										::	filename 
		integer													::	num_kpts, kpt
		real(dp)												::	raw_kpt(3), tmp
		!
		filename = 'kpts'
		inquire(file=filename, exist= read_kptsgen_pl_file)
		if( read_kptsgen_pl_file ) then
			!
			open(unit=100, file=filename, form='formatted', action='read', access='stream', status='old')
			read(100,*) num_kpts, raw_kpt(1)
			!
			allocate(	kpt_latt(3,num_kpts)	)
			!
			do kpt = 1, num_kpts
				read(100,*) raw_kpt(1:3), tmp
				kpt_latt(1:3,kpt)	= raw_kpt(1:3)
			end do
			close(100)
			write(*,'(a,i3,a,i8)')	"[#",mpi_id,"; read_kptsgen_pl_file]: success! num_kpts",num_kpts
		end if 

		!
		return 
	end function


	logical function read_geninterp_kpt_file(seed_name, kpt_latt)
		character(len=*),					intent(in)			::	seed_name
		real(dp),			allocatable,	intent(inout)		::	kpt_latt(:,:)
		character(len=50)										::	filename 
		character(len=500)										::	info_line
		integer													::	num_kpts, kpt, raw_idx
		real(dp)												::	raw_kpt(3)
		logical													::	found_file, is_frac
		!
		filename	= seed_name//'_geninterp.kpt'
		inquire(file=w90_dir//filename, exist=found_file)
		is_frac	= .false.
		!
		if(	found_file ) then
			open(unit=110, file=w90_dir//filename,	form='formatted', action='read', access='stream', 	status='old'	)
			!
			!read header
			read(110,*)
			read(110,*)	info_line
			is_frac = 	index(info_line,'frac') /=0
			
			if( is_frac ) then
				read(110,*)	num_kpts
				!
				!allocate body container
				allocate(	kpt_latt(3,	num_kpts)	)
				!
				!read body
				do kpt = 1, num_kpts
					read(110,*)		raw_idx, raw_kpt(1:3)
					kpt_latt(1:3,raw_idx)	= raw_kpt(1:3)
					if(raw_idx /= kpt)	write(*,'(a,i3,a)') '[#',mpi_id,'; read_geninterp_kpt_file]: WARNING, issues with k-mesh ordering detected'	
				end do
				write(*,'(a,i3,a,i8)')	"[#",mpi_id,"; read_geninterp_kpt_file]: success! num_kpts=",num_kpts
			else
				write(*,'(a,i3,a)')	"[#",mpi_id,"; read_geninterp_kpt_file]: file is not given in fractional coordinates (will not use file)"
			end if
			!
			close(110)
		end if
		!
		read_geninterp_kpt_file = found_file .and. is_frac
		return
	end function




	subroutine read_tb_file(seed_name, R_vect, tHopp, rHopp)
		character(len=*)									::	seed_name
		real(dp),		intent(out), allocatable			::	R_vect(:,:)
		complex(dp),	intent(out), allocatable			:: 	tHopp(:,:,:), rHopp(:,:,:,:)
		integer												::	f_nwfs, f_nSC, readLines, line, &
																cell, n
		integer												::	cell_rel(3), index(2)
		real(dp)											::	unit_cell(3,3), compl1(2),compl3(6), rTest(3), real3(3)

		

		open(unit=300, file=seed_name//'_tb.dat',form='formatted', status='old', action='read')

		!read unit cell
		read(300,*)
		read(300,*)	real3
		unit_cell(1,1:3)	= real3
		read(300,*)	 real3
		unit_cell(2,1:3) 	= real3
		read(300,*)	 real3
		unit_cell(3,1:3)	= real3
		unit_cell = unit_cell / aUtoAngstrm

		!overwrite the lattice read in from input
		a_latt	= unit_cell



		!read basis size
		read(300,*)	f_nwfs
		read(300,*)	f_nSC

		!allocate targets
		allocate(	R_vect(	3,						f_nSC		)	)
		allocate(	rHopp(	3,	f_nwfs, f_nwfs,		f_nSC		)	)
		allocate(	tHopp(		f_nwfs,	f_nwfs,		f_nSC		)	)


		
		

		!dummy read the degeneracy list
		readLines = ceiling( real(f_nSC,dp) / real(15,dp)		)
		do line = 1, readLines
			read(300,*)
		end do

		!read tHopp (eV)
		do cell = 1, f_nSC
			read(300,*)
			read(300,*) cell_rel(1:3)

			R_vect(1:3,cell)	= matmul(unit_cell,	real(cell_rel(1:3),dp)	)
			do n = 1, f_nWfs**2
				read(300,*)	index(1:2), compl1(1:2)
				!
				tHopp(index(1),index(2), cell)	= dcmplx(	compl1(1), compl1(2)		)
			end do
		end do
		

		!read rHopp (ang)
		do cell = 1, f_nSC
			read(300,*)
			read(300,*)	cell_rel(1:3)

			rTest(1:3)	= matmul(unit_cell,	real(cell_rel(1:3),dp)	)
			if(		 norm2( rTest(1:3) - R_vect(1:3,cell) )	> 1e-8_dp		) then
				write(*,*)	"[read_tb_basis]: WARNING rHopp has different sc order then tHopp"
			end if

			do n = 1, f_nWfs**2
				read(300,*)	index(1:2), compl3(1:6)

				rHopp(1,	index(1), index(2), cell)	= dcmplx( compl3(1), compl3(2)	)
				rHopp(2,	index(1), index(2), cell)	= dcmplx( compl3(3), compl3(4)	)
				rHopp(3,	index(1), index(2), cell)	= dcmplx( compl3(5), compl3(6)	)
			end do
		end do 
		!
		close(300)
		!
		!convert to a.u.
		tHopp	= 	tHopp / aUtoEv
		rHopp	= 	rHopp / aUtoAngstrm
		!
		!
		write(*,'(a,i3,a)')	"[#",mpi_id,";read_tb_file]: success!"
		return
	end subroutine 




	subroutine read_hr_file(seed_name, R_vect, H_mat)
		character(len=*),	 			intent(in) 		::	seed_name
		real(dp),		allocatable,	intent(inout)	::	R_vect(:,:)
		complex(dp),	allocatable, 	intent(inout)	::	H_mat(:,:,:)
		integer											::	f_nwfs, f_nSC, to_read, int15(15), int3(3), m,n, idx, sc, wf
		real(dp)										::	real2(2)
		integer,		allocatable						::	R_degneracy(:)

		open(unit=310, file=seed_name//'_hr.dat',	form='formatted', action='read', access='stream', status='old')

		!read header
		read(310,*)
		read(310,*)	f_nwfs
		read(310,*)	f_nSC
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
				read(310,*)		int15(1:15)
				R_degneracy(idx:idx+14)	=	int15(1:15)
				idx = idx+15	
				to_read = to_read-15
			else
				read(310,*)		int15(1:to_read)
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


			stop 'issue reading Wigner Seitz degeneracy (idx/=f_nSC) in  _hr.dat file'
		end if
		!
		!read body
		idx = 0
		do sc = 1, f_nSC
			do wf  = 1, f_nwfs**2
				!read next line
				read(310,*)		int3(1:3), m,n, 	real2(1:2)
				!
				!get Wigner Seitz vector
				if( wf==1 .and. sc==1 ) then
					R_vect(1:3,1)	= int3(1:3)
					idx = 1
				else if( .not.	is_equal_vect(	R_vect(1:3,idx),	real(int3(1:3),dp) )		) then
					idx = idx +1 
					R_vect(1:3,idx)	= int3(1:3)
					if(wf /= 1) write(*,*)	"[read_hr_file]: 	WARNING unexpected new R_vect"
				end if
				!
				!fill Hamiltonian
				H_mat(m,n, idx)	= dcmplx(	real2(1),	real2(2)	)
			end do
		end do
		close(310)
		if(idx /= f_nSC)	stop 'issue reading the body in _hr.dat file'
		!
		!
		!convert to a.u
		H_mat	= H_mat / aUtoEv
		!
		!
		
		write(*,'(a,i3,a)',advance="no")	"[#",mpi_id,";read_hr_file]: success (input interpretation: nWfs="
		write(*,'(i6,a,i6,a)')				f_nwfs, ";	nrpts=",size(R_vect,2),")!"
		return
	end subroutine



	subroutine read_r_file(seed_name, R_vect, r_mat)
		character(len=*),	 			intent(in) 		::	seed_name
		real(dp),						intent(in)		::	R_vect(:,:)
		complex(dp),	allocatable,	intent(inout)	::	r_mat(:,:,:,:)
		integer											::	f_nwfs, sc, m, n, int3(3), it
		real(dp)										::	real6(6)
		!
		open(unit=320, file=seed_name//'_r.dat',	form='formatted', action='read', access='stream', status='old')
		read(320,*)
		read(320,*)	f_nwfs
		!
		allocate(	r_mat(3, f_nwfs, f_nwfs,	size(R_vect,2)	)			)
		!
		do sc = 1, size(R_vect,2)
			do it= 1, f_nWfs**2
				read(320,*)		int3(1:3), 		m, n, 	real6(1:6)
				if( .not.	is_equal_vect(	R_vect(1:3,sc),	real(int3(1:3),dp)	)	)	then
					write(*,*)	"[read_r_file]: R_vect=",R_vect(1:3,sc)
					write(*,*)	"[read_r_file]:	input_R=",real(int3(1:3),dp)
					stop 'different R_vect order in _hr.dat and _r.dat file'
				end if
				!
				r_mat(1,m,n,sc)	=	dcmplx(	real6(1)	, real6(2)	)
				r_mat(2,m,n,sc)	=	dcmplx(	real6(3)	, real6(4)	)
				r_mat(3,m,n,sc)	=	dcmplx(	real6(5)	, real6(6)	)
			end do
		end do
		close(320)
		!
		!convert angstom to (au)
		r_mat	= r_mat	/ aUtoAngstrm
		!
		!
		return
	end subroutine





end module file_io
