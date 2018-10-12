module matrix_math

	implicit none

	private
	public						::			my_Levi_Civita, crossP,                 &
                                            zheevr_wrapper, zheevd_wrapper,         & 
                                            uni_gauge_trafo,                        &
                                            is_equal_vect,                          &
                                            is_equal_mat,                           &
                                            is_herm_mat,                            &
                                            convert_tens_to_vect,                   &
                                            blas_matmul,                            &
                                            matrix_comm                                         


    interface crossP
        module procedure   real_crossP
        module procedure   cplx_crossP
    end interface crossP

    interface  is_equal_vect
        module procedure    real_is_equal_vect
        module procedure    cplx_is_equal_vect
    end interface is_equal_vect

    interface is_equal_mat
        module procedure d_is_equal_mat
        module procedure z_is_equal_mat
    end interface is_equal_mat

    interface convert_tens_to_vect
        module procedure    real_tens_to_vect
        module procedure    cplx_tens_to_vect
    end interface convert_tens_to_vect

    interface blas_matmul
        module procedure    d_blas_matmul
        module procedure    z_blas_matmul  
    end interface blas_matmul

    interface is_herm_mat
        module procedure    z_is_herm_mat
    end interface is_herm_mat

    interface matrix_comm
        module procedure    d_matrix_comm
        module procedure    z_matrix_comm
    end interface matrix_comm








    integer,        parameter   ::  dp              = kind(0.d0)

	contains




!public:
    integer pure function my_Levi_Civita(i,j,k)
        !Hard coded Levi Civita tensor
        integer,        intent(in)      :: i,j,k
        logical                         :: even, odd
        !
        my_Levi_Civita  =   0
        !
        even    = (i==1 .and. j==2 .and. k==3) .or. (i==2 .and. j==3 .and. k==1) .or. (i==3 .and. j==1 .and. k==2)
        odd     = (i==3 .and. j==2 .and. k==1) .or. (i==1 .and. j==3 .and. k==2) .or. (i==2 .and. j==1 .and. k==3)
        !
        if(even)        then
                                my_Levi_Civita  =  1
        else if(odd)    then    
                                my_Levi_Civita  = -1
        end if
        !
        return
    end function 



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
    	call errCheck(n,info,jobz)
    	return
    end subroutine




    subroutine zheevd_wrapper(A, w)
		!solves standard (non general) eigenvalue problem with following mkl routine 
		!return eigVectors stored in A, eigValues in w
		!https://software.intel.com/en-us/node/469182        
		complex(dp) , intent(inout)   :: A(:,:)
		real(dp)    , intent(out)     :: w(:)
            !
            !set up workspace arrays (assuming jobz='V' and using zheevd)
            complex(dp), allocatable, dimension(:)   :: work
            real(dp)   , allocatable, dimension(:)   :: rwork
            integer    , allocatable, dimension(:)   :: iwork
            character(len=1)                         :: jobz,uplo
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
            !
            jobz='V'
            uplo='U' !test this behaviour
            !
            !solve eigenvalue problem
            call zheevd(jobz, uplo, n, A, n , w, work, lwork, rwork, lrwork, iwork, liwork, info)
            !
            !check if system was solved correctly
            call errCheck(n,info,jobz)
            return
	end subroutine



    subroutine uni_gauge_trafo(U_mat, M_mat)
        !   does a gauge trafo of the gauge-covariant matrix M_mat
        !   by applying the unitary matrix U_mat as shown in 
        !   Eq.(21), PRB 74, 195118, (2006)
        !
        !   since U is defined differently in the paper and in the lapack
        !   consider the following
        !
        !   want:   M^(hbar)    
        !   
        !    paper:      U^*    H^(W)   U   = H^(H)
        !           and M^(hbar) = U^* M^(W) U
        !
        !   lapack:
        !            U  H^(W)   U^*     = H^(H) 
        !   therefore here M^(hbar) = U M^(W) U^*

        complex(dp),        intent(in)      ::  U_mat(:,:)
        complex(dp),        intent(inout)   ::  M_mat(:,:)
        complex(dp),        allocatable     ::  U_dag(:,:)
        !
        allocate(U_dag(  size(U_mat,1),size(U_mat,2)  ))
        U_dag   = conjg(    transpose( U_mat )  )
        !
        !       WORKING:
        M_mat   =   blas_matmul(    M_mat,  U_mat   )
        M_mat   =   blas_matmul(    U_dag,  M_mat   )
        !
        !       NOT WORKING:   
        !M_mat   =   blas_matmul(    M_mat,  U_dag   )
        !M_mat   =   blas_matmul(    U_mat,  M_mat   )
        !       
        return
    end subroutine





!----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------
!------------INTERFACES------------------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------------------------------------
!
!
    pure function real_crossP(a,b)
        !cross product of two real 3dim vectors a,b
        real(dp), dimension(3)              :: real_crossP
        real(dp), dimension(3), intent(in)  :: a, b
        !
        real_crossP(1)   =   a(2) * b(3)     -   a(3) * b(2)  
        real_crossP(2)   =   a(3) * b(1)     -   a(1) * b(3) 
        real_crossP(3)   =   a(1) * b(2)     -   a(2) * b(1) 
        !
        return
    end function

    pure function cplx_crossP(a,b)
        !cross product of two complex 3dim vectors a,b
        complex(dp)                 :: cplx_crossP(3)
        complex(dp),    intent(in)  :: a(3), b(3)
        !
        cplx_crossP(1)   =   a(2) * b(3)     -   a(3) * b(2)  
        cplx_crossP(2)   =   a(3) * b(1)     -   a(1) * b(3) 
        cplx_crossP(3)   =   a(1) * b(2)     -   a(2) * b(1) 
        !
        return
    end function
!
!----------------------------------------------------------------------------------------------------------------------
!
!
    logical pure function real_is_equal_vect(acc, a,b)
        real(dp),     intent(in)     ::   acc
        real(dp),     intent(in)     ::   a(:), b(:)
        integer                      ::   idx
        !
        real_is_equal_vect  =   ( size(a,1) == size(b,1) ) 
        !
        if(real_is_equal_vect) then
            loop: do idx = 1, size(a,1)
                real_is_equal_vect =    real_is_equal_vect  .and.      abs(   a(idx) - b(idx)  )  < acc
                if(.not. real_is_equal_vect) exit loop
            end do loop
        end if
        !
        return
    end function

    logical pure function cplx_is_equal_vect(acc, a,b)
        real(dp),           intent(in)  ::   acc
        complex(dp),        intent(in)  ::   a(:), b(:)
        integer                         ::   idx
        !
        cplx_is_equal_vect = ( size(a,1) == size(b,1) ) 
        !
        if(cplx_is_equal_vect) then
            loop: do idx = 1, size(a,1)
                cplx_is_equal_vect =    cplx_is_equal_vect  .and.      abs(   a(idx) - b(idx)  )  < acc
                 if(.not. cplx_is_equal_vect) exit loop
            end do loop
        end if
        !
        return
    end function
!
!----------------------------------------------------------------------------------------------------------------------
!
!
    logical pure function  d_is_equal_mat(acc, A, B)
        real(dp),   intent(in)          ::  acc
        real(dp),   intent(in)          ::  A(:,:), B(:,:)
        real(dp)                        ::  delta
        integer                         ::  n, m
        !
        d_is_equal_mat  =   (   size(A,1) == size(B,1) .and. size(A,2) == size(B,2) )
        !
        if( d_is_equal_mat  ) then
            loop_out: do n = 1, size(A,2)
                do m = 1, size(A,1)
                    delta   =   abs(    A(m,n) - B(m,n) )
                    !
                    d_is_equal_mat  = d_is_equal_mat .and. (delta < acc)
                    !
                    if(.not. d_is_equal_mat) exit loop_out
                end do
            end do loop_out
        end if
        !
        return
    end function

    logical pure function z_is_equal_mat(acc, A, B)
        real(dp),       intent(in)      ::  acc
        complex(dp),    intent(in)      ::  A(:,:), B(:,:)
        real(dp)                        ::  delta
        integer                         ::  n, m
        !
        z_is_equal_mat  =   (   size(A,1) == size(B,1) .and. size(A,2) == size(B,2) )
        !
        if( z_is_equal_mat  ) then
            loop_out: do n = 1, size(A,2)
                do m = 1, size(A,1)
                    delta   =   abs(    A(m,n) - B(m,n) )
                    !
                    z_is_equal_mat  = z_is_equal_mat .and. (delta < acc)
                    !
                    if(.not. z_is_equal_mat)    exit loop_out
                end do
            end do loop_out
        end if
        !
        return
    end function
!----------------------------------------------------------------------------------------------------------------------
!
!
    logical pure function z_is_herm_mat(H)
        complex(dp),    intent(in)      ::      H(:,:)
        integer                         ::      n, m
        !
        z_is_herm_mat   =       ( size(H,1) == size(H,2)  )
        !
        if(z_is_herm_mat)   then      
            columns: do m = 1, size(H,2)
                do n = 1, size(H,1)
                    !
                    !
                    if(  n/=m       .and.       abs( H(n,m) - conjg(H(m,n)) ) >   1e-10_dp   ) then
                        z_is_herm_mat   = .false.
                        exit columns
                    end if
                    !
                    !
                end do
            end do columns 
        end if
        !
        return
    end function
!----------------------------------------------------------------------------------------------------------------------
!
!
    pure subroutine real_tens_to_vect(tens, vect)
        !   converts om_tens to a vector by applying Levi Cevita Tensor
        !   see PRB 74, 195118 (2006)   Eq.(5)
        real(dp),           intent(in)      ::  tens(3,3)
        real(dp),           intent(out)     ::  vect(3)
        integer                             ::  a, b, c
        !
        vect =  0.0_dp
        do c = 1, 3
            !
            do b = 1, 3
                do a = 1,3
                    vect(c) = vect(c) + real(my_Levi_Civita(a,b,c),dp) * tens(a,b)
                end do
            end do
            !
        end do
        return
    end subroutine

    pure subroutine cplx_tens_to_vect(tens, vect)
        !   converts om_tens to a vector by applying Levi Cevita Tensor
        !   see PRB 74, 195118 (2006)   Eq.(5)
        complex(dp),        intent(in)      ::  tens(3,3)
        complex(dp),        intent(out)     ::  vect(3)
        integer                             ::  a, b, c
        !
        vect =   cmplx(0.0_dp,0.0_dp, dp)
        do c = 1, 3
            !
            do b = 1, 3
                do a = 1,3
                    vect(c) = vect(c) + real(my_Levi_Civita(a,b,c),dp) * tens(a,b)
                end do
            end do
            !
        end do
        return
    end subroutine
!----------------------------------------------------------------------------------------------------------------------
!
!
    function d_matrix_comm(A, B) result(C)
        !   calculates the commutator of matrix A with matrix B
        !       C   = [A,B] =   AB - BA
        real(dp),    intent(in)      ::  A(:,:),     B(:,:)
        real(dp),    allocatable     ::  C(:,:)
        !
        if(size(A,2) /= size(B,1)   .or. size(B,2) /= size(A,1) ) stop "real_matrix_comm: mat sizes not matching"
        !
        !
        allocate(   C(size(A,1),size(B,2))     )
        !
        C   = blas_matmul(A, B)
        C   = C - blas_matmul(B, A)
        !
        return
    end function

    function z_matrix_comm(A,  B) result(C)
        !   calculates the commutator of matrix A with matrix B
        !       C   = [A,B] =   AB - BA
        complex(dp),    intent(in)      ::  A(:,:),     B(:,:)
        complex(dp),    allocatable     ::  C(:,:)
        !
        if(size(A,2) /= size(B,1)   .or. size(B,2) /= size(A,1) ) stop "real_matrix_comm: mat sizes not matching"
        !
        !
        allocate(   C(size(A,1),size(B,2))     )
        !
        C   = blas_matmul(A, B)
        C   = C - blas_matmul(B, A)
        !
        return
    end function
!----------------------------------------------------------------------------------------------------------------------







!private helpers
    subroutine errCheck(n,info, jobz)
        !small subroutine used in subroutine eigSolver
        implicit none
        integer    , intent(in)        :: n,info
        character(len=1), intent(in)   :: jobz
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






    function d_blas_matmul(A, B) result(C)
        real(dp),          allocatable      ::  C(:,:)
        real(dp),        intent(in)         ::  A(:,:), B(:,:)
        !
        character(len=1)                    ::  transa, transb
        integer                             ::  m, n, k, lda, ldb, ldc
        real(dp)                            ::  alpha, beta
        !
        m       =   size(A,1)
        n       =   size(B,2)
        k       =   size(A,2)
        lda     =   size(A,1)
        ldb     =   size(B,1)
        ldc     =   m
        !
        allocate(   C(m,n)     )
        !
        transa  = 'n'
        transb  = 'n'
        alpha   =  1.0_dp   
        beta    = 0.0_dp
        !
        call dgemm(transa, transb, m, n, k, alpha, A(1,1), lda, B(1,1), ldb, beta, C(1,1), ldc)
    end function




    function z_blas_matmul(A, B) result(C)
        complex(dp),        allocatable     ::  C(:,:)
        complex(dp),        intent(in)      ::  A(:,:), B(:,:)
        !
        character(len=1)                    ::  transa, transb
        integer                             ::  m, n, k, lda, ldb, ldc
        complex(dp)                         ::  alpha, beta
        !
        m       =   size(A,1)
        n       =   size(B,2)
        k       =   size(A,2)
        lda     =   size(A,1)
        ldb     =   size(B,1)
        ldc     =   m
        !
        if( k /= ldb)   stop  '[z_blas_matmul]:   warning matrix size do not match correctly'  
        !
        allocate(   C(m,n)     )
        !
        transa  = 'n'
        transb  = 'n'
        alpha   = cmplx(   1.0_dp   ,   0.0_dp  ,   dp)
        beta    = cmplx(   0.0_dp   ,   0.0_dp  ,   dp)
        !
        call zgemm(transa, transb, m, n, k, alpha,   A(1,1)    , lda,  B(1,1)    , ldb, beta,    C(1,1)      , ldc)
        !C = matmul(A,B)
    end function


end module

    



    


    

!    subroutine normalCheck(n,eigVec)
!       !check if eigVec are orthonormal this should always be the case
!       integer     , intent(in)      :: n
!       complex(dp) , intent(in)      :: eigVec(:,:)
!       integer                       :: i1,i2
!       real(dp)                      :: tmp
!       do i2=1,n
!             tmp = 0.0_dp
!             do i1=1,n
!                   tmp = tmp + dconjg(eigVec(i1,i2))*eigVec(i1,i2)
!             end do
!             if(abs(tmp)-1.0_dp>1e-12_dp) then
!                   write(*,*)"[solver/normalCheck]: WARNING eigVec not normal: ",dsqrt(tmp)
!             end if
!       end do
!       !
!       return
!   end subroutine

