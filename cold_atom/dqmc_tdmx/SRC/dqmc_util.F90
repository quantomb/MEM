module DQMC_Util
  
  use LAPACK_MOD
  use BLAS_MOD

  implicit none 
  
  ! 
  ! This module contains basic utilities used by other DQMC codes. 
  ! It also defines basic parameters.
  ! 
  ! Subroutine List
  ! ===============
  !
  !    DQMC_MatDiff(A, B) : evaluates the difference of two matrices.
  !    DQMC_Eye(A)        : returns A as an identity matrix.
  !    DQMC_ScaleCol(n, A, D, inv) : compute A*D or A*inv(D)
  !    DQMC_ScaleRow(n, A, D, inv) : compute D*A or inv(D)*A
  !    Error(message, no) : print out an error message and stop the program.
  !    ran0(n, var, seed)     : random number generators
  !    ran2(idum) result(ran) : random number generators from old program.
  !    dumpA(A, m, n, OPT): print out the content of A
  !
  ! Parameters
  ! ==========
  !
  integer,  parameter :: WP = kind(1.0d0)  ! work precision
  real(WP), parameter :: ZERO = 0.0D0      ! constant 0
  real(WP), parameter :: ONE  = 1.0D0      ! constant 1
  real(WP), parameter :: TWO  = 2.0D0      ! constant 1
  real(WP), parameter :: HALF = 0.5D0      ! constant 1

  integer,  parameter :: STDERR = 0        ! standard error output
  integer,  parameter :: STDOUT = 6        ! standard output
  integer,  parameter :: STDIN  = 5        ! standardinput

  character(*), parameter :: FMT_STRINT  = "(a30, i12)"
  character(*), parameter :: FMT_STRDBL  = "(a30, f19.6)"
  character(*), parameter :: FMT_VALERR  = "(a30, f12.6,' +- ',f12.6)"
  character(*), parameter :: FMT_INTPAR  = "(i3,i3)"
  character(*), parameter :: FMT_DBLINE  = "(76('='))"
  character(*), parameter :: FMT_SGLINE  = "(76('-'))"

contains

  !--------------------------------------------------------!
  ! Matrix computations
  !--------------------------------------------------------|

  function DQMC_MatDiff(n, A, B) result(diff)
    !
    ! Purpose
    ! =======
    !    This function computes sum(abs(A-B)).
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A and B have the same dimension.
    !    On return, A = A - B, and B is untouched.
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n         ! the order of A and B
    real(WP), intent(in)    :: A(n,n)    ! 
    real(WP), intent(in)    :: B(n,n)    !
    !
    ! ... Return value ...
    !
    real(WP) :: diff

    ! ... Blas function ...
    real(WP) :: ddot
    
    ! ... Local scalar ...
    integer  :: i

    ! ... Executable ...

    diff = ZERO
      
    do i = 1, n
       call daxpy(n, -ONE, B(1,i), 1, A(1,i), 1)
       diff = diff+ddot(n, A(1,i), 1, A(1,i), 1)
    end do

    diff = sqrt(diff)/n/n

  end function DQMC_MatDiff

  !--------------------------------------------------------!

  function DQMC_MatNorm(n, A) result(norm)
    !
    ! Purpose
    ! =======
    !    This function computes sum(abs(A-B)).
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A and B have the same dimension.
    !    On return, A = A - B, and B is untouched.
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n         ! the order of A and B
    real(WP), intent(in)    :: A(n,n)    ! 
    !
    ! ... Return value ...
    !
    real(WP) :: norm

    ! ... Local scalar ...
    integer  :: i, j
    real(WP) :: maxa, temp

    ! ... Executable ...

    norm = ZERO
    maxa = ZERO

    do i = 1, n
       do j = 1, n
          temp = abs(A(i,j))
          if (temp>maxa) then
             maxa = temp
          end if
       end do
    end do

    do i = 1, n
       do j = 1, n
          temp = A(i,j)/maxa
          norm = norm+ temp*temp
       end do
    end do

    norm = sqrt(norm)*maxa

  end function DQMC_MatNorm

  !--------------------------------------------------------!

  subroutine DQMC_Eye(n, A)
    ! 
    ! Purpose
    ! =======
    !    This subroutine returns A as an identity matrix.
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A is square.
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n        ! The order of A
    real(WP), intent(inout) :: A(n,n)   ! returned identity
   
    ! ... Local scalar ...
    integer  :: i
   
    ! ... Executable ...

    A = ZERO

    do i = 1,n
       A(i,i) = ONE
    end do
    
  end subroutine DQMC_Eye

  !--------------------------------------------------------!

  subroutine DQMC_Trans(n, At, A)
    !
    ! Purpose
    ! =======
    !    This subroutine returns the transpose of A
    !
    ! Argument
    ! ========
    integer, intent(in)     :: n
    real(wp), intent(inout) :: At(n,n)
    real(wp), intent(in)    :: A(n,n)
    
    ! ... local scalar 
    integer :: i
    
    !! decide the sgn of det(Q)
    do i = 1, n
       At(i,1:n) = A(1:n,i)
    end do
    
  end subroutine DQMC_Trans

  !--------------------------------------------------------!
  
  subroutine DQMC_ScaleCol(n, A, D)
    ! 
    ! Purpose
    ! =======
    !    This subroutine computes A = A*D.
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A is order n and D is of length n.
    !
    ! Argument
    ! ========
    integer,  intent(in)    :: n         ! The order of A
    real(WP), intent(inout) :: A(n,n)    ! 
    real(WP), intent(in)    :: D(n)      !

    ! ... Local scalar ...
    integer  :: i

    ! ... Executable ...

    ! A = A*D
    do i = 1, n
       call dscal(n, D(i), A(1, i), 1)
    end do

  end subroutine DQMC_ScaleCol
  
  !--------------------------------------------------------!
  
  subroutine DQMC_ScaleRow(n, A, D)
    ! 
    ! Purpose
    ! =======
    !    This subroutine computes A = D*A.
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A is order n and D is of length n.
    !
    ! Argument
    ! ========
    integer,  intent(in)    :: n        ! The order of A
    real(WP), intent(inout) :: A(n,n)   !
    real(WP), intent(in)    :: D(n)     ! 

    ! ... Local scalar ...
    integer  :: i
    
    ! A = D*A
    do i = 1, n
       call dscal(n, D(i), A(i,1), n)
    end do

  end subroutine DQMC_ScaleRow

  !--------------------------------------------------------!
  
  subroutine DQMC_ScaleColInv(n, A, D)
    ! 
    ! Purpose
    ! =======
    !    This subroutine computes A=A*inv(D).
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A is order n and D is of length n.
    !
    ! Argument
    ! ========
    integer,  intent(in)    :: n         ! The order of A
    real(WP), intent(inout) :: A(n,n)    ! 
    real(WP), intent(in)    :: D(n)      !

    ! ... Local scalar ...
    integer  :: i

    ! ... Executable ...

    ! A = A*inv(D)
    do i = 1, n
       call dscal(n, ONE/D(i), A(1, i), 1)
    end do

  end subroutine DQMC_ScaleColInv

  !--------------------------------------------------------!
  
  subroutine DQMC_ScaleRowInv(n, A, D)
    ! 
    ! Purpose
    ! =======
    !    This subroutine computes A = inv(D)*A.
    !
    ! Pre-assumption
    ! ==============
    !    Matrix A is order n and D is of length n.
    !
    ! Argument
    ! ========
    integer,  intent(in)    :: n        ! The order of A
    real(WP), intent(inout) :: A(n,n)   !
    real(WP), intent(in)    :: D(n)     ! 

    ! ... Local scalar ...
    integer  :: i
    
    ! ... Executable ...
    ! A = inv(D)*A
    do i = 1, n
       call dscal(n, ONE/D(i), A(i,1), n)       
    end do

  end subroutine DQMC_ScaleRowInv
 
  !--------------------------------------------------------------------!
  ! Statistics
  !--------------------------------------------------------------------!

  subroutine DQMC_SignJackKnife(n, avg, err, x, y, sgn, sum_sgn)
    !
    ! Purpose
    ! =======
    !    This subroutine implements delete-1 JackKnife method for
    !    the data and sign.
    !    X = (x1, x2, ..., xn) be the input data.
    !    sgn = (sgn1, sgn2, ..., sgnn)
    !
    !    Y = (y1, y2, ..., yn) is the Jacknife resampling of X with sign.
    !    
    !    where y_i = (sum(x)-x_i)/sgn_i
    !    The JackKnife variance of X with sign is defined as 
    !
    !     n-1 
    !    ----- sqrt(sum(y_i-avg_y)^2)
    !      n
    !
    !    where avg_y = sum(y)/n
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n
    real(wp), intent(out)   :: avg
    real(wp), intent(out)   :: err
    real(wp), intent(in)    :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in)    :: sgn(:)
    real(wp), intent(in)    :: sum_sgn
    
    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12
    
    ! ... Local variable ...
    real(wp) :: sum_x, avg_y
    
    ! ... Executable ...
    
    ! standard division
    sum_x  = sum(x)
    y      = (sum_x-x(1:n))/sgn(1:n)
    avg_y  = sum(y)/n

    y    = y - avg_y
    y    = y * y
    err  = sum(y)*(n-1)/n
    err  = sqrt(err)

    ! compute average
    avg = sum_x/sum_sgn
    
    ! If error is small enough, then regard it as 0.
    if (err .lt. TOL*abs(avg)) then
       err = ZERO
    end if
    
  end subroutine DQMC_SignJackKnife
  
  !--------------------------------------------------------------------!

  subroutine DQMC_JackKnife(n, avg, err, x, y, sgn, sum_sgn)
    !
    ! Purpose
    ! =======
    !    This subroutine implements delete-1 JackKnife method. Let
    !    X = (x1, x2, ..., xn) be the input data.
    !    Y = (y1, y2, ..., yn) is the Jacknife resampling of X.
    !    
    !    where y_i = (sum(x)-x_i)/(n-1)
    !    The JackKnife variance of X is defined as 
    !
    !     n-1 
    !    ----- sqrt(sum(y_i-avg_y)^2)
    !      n
    !
    !    where avg_y = sum(y)/n, which equals to avg_x
    !    
    ! Arguments
    ! =========
    integer,  intent(in)    :: n
    real(wp), intent(out)   :: err
    real(wp), intent(out)   :: avg
    real(wp), intent(in)    :: x(n)
    real(wp), intent(inout) :: y(n)
    real(wp), intent(inout) :: sgn(n)
    real(wp), intent(inout) :: sum_sgn


    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12
    
    ! ... Local variable ...
    real(wp) :: sum_x
    
    ! ... Executable ...

    sum_x  = sum(x)
    ! compute average
    avg    = sum_x/ n

    ! sgn and sum sgn will be used for other analysis
    sgn    = (sum_x-x(1:n))
    sum_sgn= sum_x

    ! compute y
    y    = sgn/(n-1)
    ! avg_y = avg_x 
    y    = y - avg
    y    = y * y
    err  = (sum(y)*(n-1))/n
    err  = sqrt(err)

    ! If error is small enough, then regard it as 0.
    if (err .lt. TOL*abs(avg)) then
       err = ZERO
    end if
    
  end subroutine DQMC_JackKnife
  
  !--------------------------------------------------------------------!

  subroutine DQMC_GetErr(n, err, avg, list)
    !
    ! Purpose
    ! =======
    !    This subroutine computes error of the measurements.
    !
    ! Arguments
    ! =========
    integer,  intent(in)    :: n
    real(wp), intent(out)   :: err
    real(wp), intent(in)    :: avg
    real(wp), intent(inout) :: list(n)
    
    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12
    
    ! ... Local variable ...
    real(wp) :: tmp
    
    ! ... Executable ...
    
    ! compute average
    tmp  = sum(list)/n
    
    ! standard division
    list =  list - tmp
    list = list * list
    err  = (sum(list)*(n-1))/n
    err  = sqrt(err)
    
    ! If error is small enough, then regard it as 0.
    if (err .lt. TOL*abs(avg)) then
       err = ZERO
    end if
    
  end subroutine DQMC_GetErr  

  !-----------------------------------------------------------------!

  subroutine DQMC_GetErr1(n, data, avg, err)
    !
    ! Purpose
    ! =======
    !    This subroutine computes error of the measurements.
    !
    ! Arguments
    ! =========
    integer, intent(in)   :: n
    real(wp), intent(in)  :: data(n)
    real(wp), intent(out) :: avg, err
 
    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12

    ! ... local vars ...
    integer :: i
    real(wp):: s
    
    ! Executable
    avg = sum(data) / n
    s = ZERO
    do i = 1, n
       s = s + (data(i)-avg)**2
    end do
    s = s / n
    err = sqrt(s)/sqrt(n-ONE)

    if (err .lt. TOL*abs(avg)) then
       err = ZERO
    end if

  end subroutine DQMC_GetErr1
   
  !--------------------------------------------------------!

  subroutine DQMC_GetErr2(n, sm, ssq, avg, err)
    !
    ! Purpose
    ! =======
    !    This subroutine computes error of the measurements.
    !
    ! Arguments
    ! =========
    integer, intent(in)   :: n
    real(wp), intent(in)  :: sm, ssq
    real(wp), intent(out) :: avg, err
 
    ! ... Parameter ...
    real(wp), parameter :: TOL = 1.0D-12

    ! ... local vars ...
    real(wp):: s
    
    ! Executable
    avg = sm / n
    s   = ssq / n - avg*avg
    err = sqrt(s) / sqrt(n-ONE)

    if (err .lt. TOL*abs(avg)) then
       err = ZERO
    end if

  end subroutine DQMC_GetErr2

  !--------------------------------------------------------!
  ! Miscellanceous
  !--------------------------------------------------------!

  subroutine DQMC_Error(message, no)
    ! 
    ! Purpose
    ! =======
    !    This subroutine prints out an error message and
    !    a message number, and then stop the program.
    !
    ! Arguments
    ! =========
    ! 
    character(*), intent(in) :: message   ! Error message
    integer, intent(in)      :: no        ! Message number

    ! ... Executable ...

    write(STDERR,*) "Error: ", message
    stop

  end subroutine DQMC_Error
  
  !--------------------------------------------------------!

  subroutine DQMC_Warning(message, no)
    ! 
    ! Purpose
    ! =======
    !    This subroutine prints out an error message and
    !    a message number.
    !
    ! Arguments
    ! =========
    ! 
    character(*), intent(in) :: message   ! Warning message
    integer, intent(in)      :: no        ! Message number

    ! ... Executable ...

    write(STDERR,*) "Warning: ", message

  end subroutine DQMC_Warning
  
  !--------------------------------------------------------!

  subroutine ran0(n, var, seed)
    ! 
    ! Purpose
    ! =======
    !    Random number generator. This subroutine uses 
    !    LAPACK's random number generator dlaruv or 
    !    SPRNG to generate a list of random numbers
    !
    ! SPRNG DECLARE
    ! =============
#define SIMPLE_SPRNG
#define FLOAT_GEN
#include "sprng_f.h" 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)     :: n          ! length of the list
    real(wp), intent(out)   :: var(n)     ! random number to return
    integer, intent(inout)  :: seed(4)    ! random seeds

    ! ... local scalar ...
    integer, parameter :: max_len = 128   ! This max length is 
                                          ! defined by dlaruv
    integer :: i
    ! ... Executable ...

    var(1:n) = ZERO

    ! USE sprng instead of lapack code dlaruv
    do i = 1, n
       var(i) = sprng()
    end do

!    call dlarnv(DLARNV_UNI_0_1, seed, n, var)

  end subroutine ran0

  !--------------------------------------------------------!

  subroutine ran1(n, var, seed)
    ! 
    ! Purpose
    ! =======
    !    Random number generator. This subroutine uses 
    !    LAPACK's random number generator dlaruv or 
    !    SPRNG to generate a list of random numbers
    !
    ! SPRNG DECLARE
    ! =============
#define SIMPLE_SPRNG
#define FLOAT_GEN
#include "sprng_f.h"
    !
    ! Arguments
    ! =========
    ! 
    integer, intent(in)     :: n          ! length of the list
    real(wp), intent(out)   :: var(n)     ! random number to return
    integer, intent(inout)  :: seed(4)    ! random seeds

    ! ... local scalar ...
    integer, parameter :: max_len = 128   ! This max length is 
                                          ! defined by dlaruv
    integer :: i
    ! ... Executable ...

    var(1:n) = ZERO
    ! USE SPRNG
    do i = 1, n
       var(i) = TWO*sprng()-ONE
    end do

!    call dlarnv(2, seed, n, var)

  end subroutine ran1

 !--------------------------------------------------------!

!  subroutine ranN(n, var, seed)
    ! 
    ! Purpose
    ! =======
    !    Random number generator. This subroutine uses 
    !    LAPACK's random number generator dlaruv to generate
    !    a list of random numbers
    !
    ! Arguments
    ! =========
    ! 
!    integer, intent(in)     :: n          ! length of the list
!    real(wp), intent(out)   :: var(n)     ! random number to return
!    integer, intent(inout)  :: seed(4)    ! random seeds

!    ! ... local scalar ...
!    integer, parameter :: max_len = 128   ! This max length is 
                                          ! defined by dlaruv
!    ! ... Executable ...

!    var(1:n) = ZERO
!    call dlarnv(3, seed, n, var)

!  end subroutine ranN

  !-------------------------------------------------------------!

  subroutine dumpA(m, n, A, OPT)
    implicit none
    ! 
    ! Purpose
    ! =======
    !    This subroutine prints the content of matrix A.
    !    *** for internal debugging only
    !
    ! Arguments
    ! =========
    ! 
    integer, intent(in)  :: m, n          ! dimension of A
    real(wp), intent(in) :: A(1:m,1:n)    ! matrix A
    integer, intent(in)  :: OPT           ! output device

    ! ... Local variables ...
    character(20) fmt
    integer i 
    real(wp) :: temp(n)

    ! ... Executable ...

    write(fmt,"(A,I3,A)") "(",n,"F20.15)"
    
    do i=1, m
       temp = A(i,1:n)
       write (OPT, fmt) temp
    end do
    
  end subroutine dumpA

  !--------------------------------------------------------!
    
  subroutine DQMC_Print_RealArray(n, m, title, label, avg, err, OPT)
    !
    ! Purpose
    ! =======
    !    This subroutine prints out the content of avg +- err
    !   
    ! Arguments
    ! =========
    integer,  intent(in)    :: OPT, n, m
    character(*), intent(in):: title, label(:)
    real(wp), intent(in)    :: avg(:,:), err(:,:)
    
    ! ... Local variable ...
    integer  :: i, j
    
    ! ... Executable ...
    
    write(OPT,*) title
    if (n .gt. 0) then
       do i = 1, n
          write(OPT,*) label(i)
          do j = 1, m
             write(OPT, "(i3,f10.6,' +-',f10.6)") j-1, avg(i,j), err(i,j)
          end do
       end do
    else
       do j = 1, m
          write(OPT, "(a30,f10.6,' +-',f10.6)") label(j), avg(j,1), err(j,1)
       end do
    end if
    write(OPT,FMT_DBLINE)

  end subroutine DQMC_Print_RealArray

  !--------------------------------------------------------------------!
  
  subroutine DQMC_Print_ComplexArray(n, m, title, label, avg, err, OPT)
    !
    ! Purpose
    ! =======
    !    This subroutine prints out the content of avg +- err
    !   
    ! Arguments
    ! =========
    integer,  intent(in)    :: OPT, n, m
    character(*), intent(in):: title, label(:)
    real(wp), intent(in)    :: avg(:,:), err(:,:)
    
    ! ... Local variable ...
    integer  :: i, j
    
    ! ... Executable ...
    
    write(OPT,*) title

    if (n .gt. 0) then
       do i = 1, n
          write(OPT,*) label(i)
          do j = 1, m
             write(OPT,"(i3,'( ',f10.6,' +-',f10.6,' ) &
                  & +i( ',f10.6,' +-',f10.6,' )')")&
                  j-1, avg(i,j), err(i,j), avg(i,j+m), err(i,j+m)
          end do
       end do
    else
       do j = 1, m
          write(OPT,"(a30,'( ',f10.6,' +-',f10.6,' ) &
               & +i( ',f10.6,' +-',f10.6,' )')") &
               label(j), avg(j,1), err(j,1), avg(j+m,1), err(j+m,1)
       end do
    end if
    
    write(OPT,FMT_DBLINE)

  end subroutine DQMC_Print_ComplexArray  

  !--------------------------------------------------------------------!

end module DQMC_Util
