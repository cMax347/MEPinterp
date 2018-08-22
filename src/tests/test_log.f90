module test_log

	implicit none


	private
	public						::		init_outFile, 		&
										write_test_results,	&
										push_to_outFile



	character(len=2)			::		out_dir="./"
	character(len=30)			::		fpath	
	integer						::		unit = 111
	logical						::		initialized=.false., print_cli



	contains




	subroutine	init_outFile(cli_print, fname)
		logical,				intent(in)				::	cli_print
		character(len=*),		intent(in)				::	fname
		character(8)  :: date
    	character(10) :: time
    	character(5)  :: zone
    	integer,dimension(8) :: values
    	!
    	print_cli	= cli_print
    	!set fpath
		fpath	=	out_dir//fname
		!
		!print header to file
		open(unit=unit, file=fpath, form='formatted', 	action='write', access='stream',	status='replace')
		call date_and_time(date,time,zone,values)
		write(unit,*)	"# MEPinterp TESTS performed on ",date," ", time," ", zone
		write(unit,*)
		close(unit)
		!
		initialized	= .true.
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
		!
		if(size(passed) /= size(label)	)	stop "non matching passed and label array"
		nTests		= size(passed)
		succ_cnt	= 0
		!
		if( .not. initialized	)	call init_outFile(.false.,"test.log")
		open(unit=unit,file=fpath, form="formatted",  action='write', status='old', position='append')
		!
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
		!
		if( .not. initialized	)	call init_outFile(.false.,"test.log")
		!
		open(unit=unit,file=fpath, form="formatted",  action='write', status='old', position='append')
			write(unit,*)	push_notification
		close(unit)
		!
		if(	print_cli	)	write(*,*)	push_notification
		!
		return
	end subroutine








end module 