module test_log

	implicit none


	private
	public						::		init_outFile, push_to_outFile



	character(len=2)			::		out_dir="./"
	character(len=30)			::		fpath	
	integer						::		unit = 111
	logical						::		initialized=.false.



	contains




	subroutine	init_outFile(fname)
		character(len=*),		intent(in)				::	fname
		character(8)  :: date
    	character(10) :: time
    	character(5)  :: zone
    	integer,dimension(8) :: values
    	!
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


	subroutine push_to_outFile(push_notification)
		character(len=*),	intent(in)		::	push_notification
		!
		if( .not. initialized	)	call init_outFile("test.log")
		!
		open(unit=unit,file=fpath, form="formatted",  action='write', status='old', position='append')
		write(unit,*)	push_notification
		close(unit)
		!
		return
	end subroutine








end module 