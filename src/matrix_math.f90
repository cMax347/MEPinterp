module matrix_math



	implicit none



	private
	public						::			zheevr_wrapper, zheevd_wrapper,       & 
                                            matrix_comm,                          &
                                            crossP, is_equal_vect


    interface crossP
        module procedure crossPreal
        module procedure crossPcplx
    end interface crossP


    integer,        parameter   ::  dp              = kind(0.d0)




	contains




!public:
	subroutine zheevr_wrapper(a, w ,z, m)
    	!https://software.intel.com/en-us/mkl-developer-reference-fortran-heevr#6ADF761A-127A-4C9B-9A2A-1A8AA4602CE1
    	!with ability to solve only selected eigenvalues (& vectors)
    	complex(dp),	intent(inout)			:: 	a(:,:)				!Hermitian matrix to be solved
    	real(dp),		intent(out)				:: 	w(:)				!selected eigenvalues in ascending order
    	complex(dp),	intent(out)				:: 	z(:,:)				!each column represents an eigenvector (	columns 1,2,...,m	)
    	integer,		intent(out)				:: 	m					!number of eigenvalues found
    	character(len=1)	 					:: 	jobz, range, uplo
    	integer									:: 	n, lda, il, iu, ldz, lwork, lrwork, liwork,  info 
    	integer,		allocatable				:: 	isuppz(:) , iwork(:)
    	real(dp)								:: 	vl, vu, abstol
    	real(dp),		allocatable				:: 	rwork(:)
    	complex(dp),	allocatable				:: 	work(:)
    	!
    	jobz	= 'V'			!also compute eigenvectors
    	range	= 'I'
    	uplo	= 'U'
    	n 		= size(a,1)
    	lda		= n
    	vl		= 0.0_dp
    	vu  	= 0.0_dp
    	il		= 1
    	iu		= size(z,2)		!number of eigenvalues to compute
    	abstol	= 1e-15_dp
    	ldz		= n
    	lwork	= 2*n
    	lrwork	= 24*n
    	liwork	= 10*n
    	if( size(z,1)/= ldz ) write(*,*)"[eigSolver2]: z array has wrong size"
    	if( size(w)/= n) write(*,*)"[eigSolver2]; w array has wrong size"
    	!
    	allocate( isuppz(2*iu)	)
    	allocate(  work(lwork)	)
    	allocate( rwork(lrwork)	)
    	allocate( iwork(liwork)	)	
    	
    	!
    	call zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &
    					m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
    	!
    	!
    	call errCheck(n,info,jobz)
    	!
    	!
    	return
    end subroutine




    subroutine zheevd_wrapper(A, w)
		!solves standard (non general) eigenvalue problem with following mkl routine 
		!return eigVectors stored in A, eigValues in w
		!https://software.intel.com/en-us/node/469182        
		complex(dp) , intent(inout)   :: A(:,:)
		real(dp)    , intent(out)     :: w(:)
            

            !set up workspace arrays (assuming jobz='V' and using zheevd)
            complex(dp), allocatable, dimension(:)   :: work
            real(dp)   , allocatable, dimension(:)   :: rwork
            integer    , allocatable, dimension(:)   :: iwork
            character*1                              :: jobz,uplo
            integer                                  :: n, info,lwork,lrwork,liwork
            n		= size(A,1)
            if(n /= size(A,2)) then
            	write(*,*)"[eigSolver]: WARNING the matrix to solve is not a square matrix"
            end if
            lwork  	=   n*n + 2*n
            lrwork 	= 2*n*n + 5*n + 1
            liwork 	=         5*n + 3
            allocate(  work( lwork) )
            allocate( rwork(lrwork) )
            allocate( iwork(liwork) )

            jobz='V'
            uplo='U' !test this behaviour

            !solve eigenvalue problem
            !call zheevd(jobz, uplo, n, a,lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
             call zheevd(jobz, uplo, n, A, n , w, work, lwork, rwork, lrwork, iwork, liwork, info)
            
            !check if system was solved correctly
            call errCheck(n,info,jobz)
            !
            if(.true.) then
                  call normalCheck(n,A)
            end if
            !
            return
	end subroutine




    subroutine  matrix_comm(A,  B,  C)
        !   calculates the commutator of matrix A with matrix B
        !       C   = [A,B] =   AB - BA
        !
        complex(dp),    intent(in)      ::  A(:,:),     B(:,:)
        complex(dp),    intent(out)     ::  C(:,:)
        !
        C   = matmul(A, B)
        C   = C - matmul(B, A)
        !
        return
    end subroutine




!private:

	subroutine errCheck(n,info, jobz)
		!small subroutine used in subroutine eigSolver
		implicit none
		integer    , intent(in) :: n,info
		character*1, intent(in) :: jobz
		!
		if(info > 0) then
		      write(*,*) '[solver/errCheck]: WARNING, Problem solving the eigenvalue problem: '
		      if(jobz .EQ. 'N') write(*,*) '[solver/errCheck]: the algorithm failed to converge; ', info ,&
		      				' off-diagonal elements of an intermediate tridiagonal form did not converge to zero;'
		      if(jobz .EQ. 'V') write(*,*) '[solver/errCheck]: the algorithm failed to compute an',&
		    ' eigenvalue while working on the submatrix lying in rows and columns ', info/(n+1), ' through ', mod(info, n+1)
		elseif(info < 0) then
		     write(*,*) '[solver/errCheck]: the ',info,'-th parameter had an illegal value.' 
		endif
		!
		return
    end subroutine


    subroutine normalCheck(n,eigVec)
    	!check if eigVec are orthonormal this should always be the case
		integer     , intent(in)      :: n
		complex(dp) , intent(in)      :: eigVec(:,:)
		integer                       :: i1,i2
		real(dp)                      :: tmp
		do i2=1,n
		      tmp = 0.0_dp
		      do i1=1,n
		            tmp = tmp + dconjg(eigVec(i1,i2))*eigVec(i1,i2)
		      end do
		      if(abs(tmp)-1.0_dp>1e-12_dp) then
		            write(*,*)"[solver/normalCheck]: WARNING eigVec not normal: ",dsqrt(tmp)
		      end if
		end do
		!
		return
	end subroutine

    logical function is_equal_vect(acc, a,b)
        real(dp),                   intent(in)     ::   acc
        real(dp),   dimension(3),   intent(in)     ::   a, b
        integer                                    ::   x
        !
        is_equal_vect = .true.
        do x = 1, 3
            is_equal_vect =    is_equal_vect  .and.      abs(   a(x) - b(x)  )  < acc
        end do
        !
        return
    end function

    function crossPreal(a,b)
        !cross product of two real 3dim vectors a,b
        real(dp), dimension(3)              :: crossPreal
        real(dp), dimension(3), intent(in)  :: a, b
        !
        crossPreal(1)   =   a(2) * b(3)     -   a(3) * b(2)  
        crossPreal(2)   =   a(3) * b(1)     -   a(1) * b(3) 
        crossPreal(3)   =   a(1) * b(2)     -   a(2) * b(1) 
        !
        return
    end function


    function crossPcplx(a,b)
        !cross product of two complex 3dim vectors a,b
        complex(dp),    intent(in)  :: a(3), b(3)
        complex(dp)                 :: crossPcplx(3)
        !
        crossPcplx(1)   =   a(2) * b(3)     -   a(3) * b(2)  
        crossPcplx(2)   =   a(3) * b(1)     -   a(1) * b(3) 
        crossPcplx(3)   =   a(1) * b(2)     -   a(2) * b(1) 
        !
        return
    end function




end module