module matrix_math

	implicit none

	private
	public						::	    get_levi_civita,                        &
                                        crossP,                                 &
                                        zheevr_wrapper,                         &
                                        zheevd_wrapper,                         &
                                        zheevx_wrapper,                         &
                                        is_equal_vect,                          &
                                        is_equal_mat,                           &
                                        is_herm_mat,                            &
                                        is_skew_herm_mat,                       &
                                        convert_tens_to_vect,                   &
                                        blas_matmul,                            &
                                        matrix_comm,                            &
                                        get_linspace                                        


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

    interface is_skew_herm_mat
        module procedure    z_is_skew_herm_mat
    end interface is_skew_herm_mat 

    interface matrix_comm
        module procedure    d_matrix_comm
        module procedure    z_matrix_comm
    end interface matrix_comm




    save



    integer,        parameter   ::  dp              = kind(0.d0)


	contains




!public:
    pure subroutine get_linspace(min, max, n_mesh, linspace)
        real(dp),                      intent(in)      ::   min, max
        integer,                       intent(in)      ::   n_mesh
        real(dp),   allocatable,       intent(out)     ::   linspace(:)
        real(dp)                                       ::   delta
        integer                                        ::   idx 
        !
        !
        if( n_mesh > 1) then
            !   
            !   RETURN LIST
            allocate(   linspace(n_mesh)  )
            delta   =   (   max -   min)    /   real(n_mesh-1,dp)    
            !
            do idx  =   1, n_mesh
                linspace(idx)   =       min     +  real(idx-1,dp)    * delta    
            end do
        else
            !
            !   RETURN SCALAR
            allocate(   linspace(1) )  
            linspace    =   min
        end if
        !
        !
        return
    end subroutine



    integer pure function my_Levi_Civita(i,j,k)
        !Hard coded Levi Civita tensor
        integer,        intent(in)      :: i,j,k
        logical                         :: even, odd
        !
        !
        even    = (i==1 .and. j==2 .and. k==3) .or. (i==2 .and. j==3 .and. k==1) .or. (i==3 .and. j==1 .and. k==2)
        odd     = (i==3 .and. j==2 .and. k==1) .or. (i==1 .and. j==3 .and. k==2) .or. (i==2 .and. j==1 .and. k==3)
        !
        if(even)        then
                                my_Levi_Civita  =  1
        else if(odd)    then    
                                my_Levi_Civita  = -1
        else 
                                my_Levi_Civita  =  0
        end if
        !
        return
    end function 


    pure subroutine get_levi_civita(lc_tens)
        integer,    intent(out)   ::  lc_tens(3,3,3)
        !
        lc_tens         =    0
        !
        lc_tens(1,2,3)  =   +1
        lc_tens(2,3,1)  =   +1
        lc_tens(3,1,2)  =   +1
        !
        lc_tens(1,3,2) =   -1
        lc_tens(3,2,1) =   -1
        lc_tens(2,1,3) =   -1        
        !
        return
    end subroutine





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
        if( info /= 0)  stop 'zheevr'
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
            if( info /= 0)  stop 'zheevd'
            return
	end subroutine




    subroutine zheevx_wrapper(A, w)
        complex(dp),        intent(inout)       ::  A(:,:)
        real(dp),           intent(in)          ::  w(:)
        real(dp)                                ::  vl, vu
        integer                                 ::  il, iu, lwork, lrwork, liwork, info, num_wann
        complex(dp),        allocatable         ::  z(:,:), work(:) 
        real(dp),           allocatable         ::  rwork(:)
        integer,            allocatable         ::  iwork(:), ifail(:)       
        !
        if (size(A,1)== size(A,2)) then
            num_wann    =   size(A,1)
        else
            stop "zheevx_wrapper was given non symmetric matrix!"
        end if
        !
        lwork=12*num_wann
        lrwork=17*num_wann
        liwork=15*num_wann
        allocate(   ifail(          num_wann            ))
        allocate(   z(             num_wann, num_wann   ))
        allocate( rwork(lrwork),    work(lwork),    iwork(liwork) )
        !
        !
        call zheevx('V', 'A', 'U', size(A,1), A(:,:), size(A,1), vl, vu, il, iu, 0.0_dp, size(w,1),& 
                            w(:), z(:,:), size(A,1), work, lwork, rwork, iwork, ifail, info)
        !
        !   return eigenvectors
        A(:,:)  = z(:,:)
        !
        !   debug
        call errCheck(num_wann,info,'V')
        if( info /= 0)  stop 'zheevx'
        return
    end subroutine


!    subroutine uni_gauge_trafo(U_mat, M_mat)
!        !   does a gauge trafo of the gauge-covariant matrix M_mat
!        !   by applying the unitary matrix U_mat as shown in 
!        !   Eq.(21), PRB 74, 195118, (2006)
!        !
!        !   since U is defined differently in the paper and in the lapack
!        !   consider the following
!        !
!        !   want:   M^(hbar)    
!        !   
!        !    paper:      U^*    H^(W)   U   = H^(H)
!        !           and M^(hbar) = U^* M^(W) U
!        !
!        !   lapack:
!        !            U  H^(W)   U^*     = H^(H) 
!        !   therefore here M^(hbar) = U M^(W) U^*
!
!        complex(dp),        intent(in)      ::  U_mat(:,:)
!        complex(dp),        intent(inout)   ::  M_mat(:,:)
!        complex(dp),        allocatable     ::  U_dag(:,:)
!        !
!        allocate(U_dag(  size(U_mat,1),size(U_mat,2)  ))
!        U_dag   = conjg(    transpose( U_mat )  )
!        !
!        !       WORKING:
!       ! M_mat   =   blas_matmul(    M_mat,  U_mat   )
!       ! M_mat   =   blas_matmul(    U_dag,  M_mat   )
!        !
!        M_mat   =   blas_matmul(    blas_matmul(U_dag, M_mat),  U_mat)
!        !       NOT WORKING:   
!        !M_mat   =   blas_matmul(    M_mat,  U_dag   )
!        !M_mat   =   blas_matmul(    U_mat,  M_mat   )
!        !       
!        return
!    end subroutine





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
        real(8),       intent(in)      ::  acc
        complex(8),    intent(in)      ::  A(:,:), B(:,:)
        real(8)                        ::  delta
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
    logical function z_is_herm_mat(H, max_err)
        complex(dp),    intent(in)      ::      H(:,:)
        real(dp),       intent(out)     ::      max_err
        real(dp)                        ::      dH_nm
        integer                         ::      n, m
        !
        max_err         =       -1.0_dp
        z_is_herm_mat   =       ( size(H,1) == size(H,2)  )
        !
        if(z_is_herm_mat)   then      
            columns: do m = 1, size(H,2)
                do n = 1, size(H,1)
                    !
                    dH_nm   =  abs( H(n,m) - conjg(H(m,n)) )  
                    !
                    if(  n/=m       .and.  dH_nm     >   1e-6_dp   ) then
                        z_is_herm_mat   = .false.
                        !write(*,*)  '[z_is_herm_mat]: herm. criterion violated by ', abs( H(n,m) - conjg(H(m,n)) )
                        !exit columns
                        if(dH_nm > max_err)     max_err =   dH_nm
                    end if
                    !
                    !
                end do
            end do columns 
        else
            write(*,*)  '[z_is_herm_mat]: WARNING  matrix is non symmetric...'
        end if
        !
        return
    end function
!----------------------------------------------------------------------------------------------------------------------
!
!
    logical function z_is_skew_herm_mat(H, max_err)
        complex(dp),    intent(in)      ::      H(:,:)
        real(dp),       intent(out)     ::      max_err
        real(dp)                        ::      dH_nm
        integer                         ::      n, m
        !
        max_err              =      -1.0_dp
        z_is_skew_herm_mat   =      ( size(H,1) == size(H,2)  )
        !
        if(z_is_skew_herm_mat)   then      
            columns: do m = 1, size(H,2)
                do n = 1, size(H,1)
                    !
                    dH_nm   =   abs( H(n,m) + conjg(H(m,n)) )
                    !
                    if(  n/=m       .and.   dH_nm     >   1e-6_dp   ) then
                        z_is_skew_herm_mat   = .false.
                        !exit columns
                        if(dH_nm > max_err)     max_err =   dH_nm
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
        real(dp),    allocatable     ::  C(:,:),     tmp(:,:)
        !
        if(size(A,2) /= size(B,1)   .or. size(B,2) /= size(A,1) ) stop "real_matrix_comm: mat sizes not matching"
        !
        !
        allocate(   C(  size(A,1),size(B,2))     )
        allocate(   tmp(size(A,1),size(B,2))     )

        !
        call blas_matmul(A, B,      C   )
        call blas_matmul(B, A,      tmp )
        C   = C - tmp
        !
        return
    end function

    function z_matrix_comm(A,  B) result(C)
        !   calculates the commutator of matrix A with matrix B
        !       C   = [A,B] =   AB - BA
        complex(dp),    intent(in)      ::  A(:,:),     B(:,:)
        complex(dp),    allocatable     ::  C(:,:),     tmp(:,:)
        !
        if(size(A,2) /= size(B,1)   .or. size(B,2) /= size(A,1) ) stop "real_matrix_comm: mat sizes not matching"
        !
        !
        allocate(   C(size(A,1),size(B,2))     )
        allocate(   tmp(size(A,1),size(B,2))   )
        !
        call blas_matmul(A, B,      C   )
        call blas_matmul(B, A,      tmp )
        C   = C - tmp
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






    subroutine d_blas_matmul(A, B, C) 
        !   C := alpha*op(A)*op(B) + beta*C
        real(dp),        intent(in)         ::  A(:,:), B(:,:)
        real(dp),        intent(out)        ::  C(:,:)
        integer                             ::  m, n, k
        !
        m       =   size(A,1)
        n       =   size(B,2)
        k       =   size(A,2)
        if( k /= size(B,1))   stop  '[d_blas_matmul]:   warning matrix size do not match correctly'  

        !
        !allocate(   C(m,n)     )
        C       =   0.0_dp
        !
       
        !
        !        CALL DGEMM('N','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)
        call dgemm('N', 'N', m, n, k, 1._dp, A(:,:), m, B(:,:), k, 0._dp, C(:,:), m)
        !C = matmul(A,B)
        return
    end subroutine




    subroutine z_blas_matmul(A, B, C) 
        !   C := alpha*op(A)*op(B) + beta*C
        complex(dp),        intent(in)      ::  A(:,:), B(:,:)
        complex(dp),        intent(out)     ::  C(:,:)
        integer                             ::  m, n, k
        !
        m       =   size(A,1)
        k       =   size(A,2)
        n       =   size(B,2)
        if( k /= size(B,1))   stop  '[z_blas_matmul]:   warning matrix size do not match correctly'  
        !
        !
        C   =   cmplx(0.0_dp,0.0_dp,dp)
        !
        call zgemm('N', 'N', m, n, k, 1._dp, A(:,:), m,  B(:,:), k, 0,    C(:,:), m)
        !C = matmul(A,B)
        !
        return
    end subroutine


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

