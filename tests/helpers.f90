!
!	module providing usefull routines for testing


module helpers	
	use	constants,			only:		dp, fp_acc

	implicit none

	private	
	public						::		&
										my_exit,			&
										!	LOG FILE
										init_outFile, 		&
										write_test_results,	&
										push_to_outFile,	&
										!	MATRIX
										random_matrix,		&
										random_vector


	character(len=30)			::		fpath="./test.log"	
	integer						::		unit = 111




	interface random_matrix
		module procedure ds_rand_mat
		module procedure zh_rand_mat	
	end interface random_matrix

	interface random_vector
		module procedure d_rand_vec
		module procedure z_rand_vec
	end interface random_vector






contains
	subroutine	my_exit(normal_exit)
		logical, 	intent(in)		::	normal_exit
		!
		if( normal_exit)	call exit(0)
		call exit(1)
		!
	end subroutine




	!
	!
	!-----------------------------------------------------------------------------------------------------------------------|
	! 			IO - HELPERS																								|
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
	subroutine	init_outFile()
		character(8)  :: date
    	character(10) :: time
    	character(5)  :: zone
    	integer,dimension(8) :: values
    	
		!print header to file
		open(unit=unit, file=fpath, form='formatted', 	action='write', status='replace')
		call date_and_time(date,time,zone,values)
		write(unit,*)	'*'
		write(unit,*)	'*'
		write(unit,*)	'*'
		write(unit,*)	'*'
		write(unit,*)	'*'
		write(unit,*)	'*'
		write(unit,*)	"# MEPinterp TESTS performed on ",date," ", time," ", zone
		write(unit,*)
		close(unit)
		!
		return
	end subroutine

	subroutine write_test_results(descriptor,	passed, label, succ_cnt, nTests)
		character(len=*),	intent(in)			::	descriptor, label(:)
		logical,			intent(in)			::	passed(:)
		integer,			intent(out)			::	succ_cnt, nTests
		integer									::	test
		character(len=10)						:: 	result
		character(len=121)						::	succ_message
		logical									::	initialized
		!
		if(size(passed) /= size(label)	)	stop "non matching passed and label array"
		nTests		= size(passed)
		succ_cnt	= 0
		!
		inquire(file=fpath, exist= initialized)
		if( .not. initialized	)	call init_outFile()
		!
		!
		open(unit=unit,file=fpath, form="formatted",  action='write', status='old', position='append')
		do test = 1, nTests
			if( passed(test) )	then 
				result	=	"passed"
				succ_cnt=	1 + succ_cnt
			else
				result	=	"not passed"
			end if
			write(succ_message,'(a,a,a,a,a,a,a)')	'[',descriptor,']: ',result,' test of "',trim(label(test)),'"'
			!
			!WRITE TO LOG FILE
			write(unit,*)	succ_message
			!call push_to_outFile(succ_message)
		end do
		write(succ_message,'(a,a,a,i1,a,i1,a)')		'[',descriptor,']: passed ',succ_cnt,' of ',nTests,' tests'
		!call push_to_outFile(succ_message)
		write(unit,*)		succ_message
		!
		return
	end subroutine

	subroutine push_to_outFile(push_notification)
		character(len=*),	intent(in)		::	push_notification
		logical								::	initialized
		!
		inquire(file=fpath, exist= initialized)
		if( .not. initialized	)	call init_outFile()
		!
		open(unit=unit,file=fpath, form="formatted",  action='write', status='old', position='append')
			write(unit,*)	trim(push_notification)
		close(unit)
		!
		return
	end subroutine



	!
	!
	!-----------------------------------------------------------------------------------------------------------------------|
	! 			MATRIX HELPERS																								|
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

	subroutine ds_rand_mat(Mat)
		real(dp), intent(out)	::	Mat(:,:)
		integer						::	n, m
		real(dp)					::	re_rand
		!
		!
		do n = 1, size(Mat,2)
			do m = 1, size(Mat,1)
				call random_number(re_rand)	
				!
				re_rand	= (re_rand - .5_dp) 							
				!
				Mat(n,m)	=	re_rand
				if(n<m) then
					Mat(m,n)	=	re_rand
				end if
				!
				!
			end do
		end do
		!
		return
	end subroutine	

	subroutine zh_rand_mat(Mat)
		complex(dp), intent(out)	::	Mat(:,:)
		integer						::	n, m
		complex(dp)					::	rand
		real(dp)					::	re_rand, im_rand
		!
		do n = 1, size(Mat,2)
			do m = 1, size(Mat,1)
				call random_number(re_rand)
				re_rand	= (re_rand - .5_dp) 		
				!
				if(	n==m ) then
					Mat(n,n)	=	re_rand
				else
					call random_number(im_rand)	
					im_rand	= (im_rand - .5_dp)					
					rand	= cmplx(	re_rand, im_rand,	dp)
					!
					Mat(n,m)	=	rand
					Mat(m,n)	= conjg(rand)
				end if
				!
				!
			end do
		end do
		!
		return
	end subroutine	
!------------------------------------------------------------------------------------------------------

	subroutine d_rand_vec(v)
		real(dp),	intent(out)		::	v(:)
		integer						:: idx
		real(dp)					:: rand, scal
		!
		call random_number(scal)
		scal	=	scal * size(v,1)
		!
		do idx = 1, size(v,1)
			call random_number(rand)
			rand = (rand-.5_dp)	* scal
			!
			!
			v(idx)	= rand
		end do
		!
		return
	end subroutine

	subroutine z_rand_vec(v)
		complex(dp),	intent(out)	::	v(:)
		integer						:: 	idx
		real(dp)					:: 	rand_re, rand_im, scal
		complex(dp)					::	rand
		!
		call random_number(scal)
		scal	=	scal * size(v,1)
		!
		do idx = 1, size(v,1)
			call random_number(rand_re)
			call random_number(rand_im)
			!
			rand_re = 	(rand_re-.5_dp)	* scal
			rand_im = 	(rand_im-.5_dp)	* scal
			rand 	=	cmplx(	rand_re	,	rand_im	,	dp)
			!
			!
			v(idx)	= rand
		end do
		!
		return
	end subroutine




end module