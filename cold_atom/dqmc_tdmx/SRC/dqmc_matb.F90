module DQMC_MATB
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_Struct
  use BLAS_MOD
  use LAPACK_MOD
  use DQMC_WSPACE

  implicit none 

  ! 
  ! This module defines the type and subroutines for propagator B.
  !
  !  Data Type
  !  =========
  !
  type matB
     integer  :: n                      ! dim of B
     real(wp), pointer :: B(:,:)        ! Matrix B 
     real(wp), pointer :: Bi(:,:)       ! Inverse of matrix B
     character(12)     :: name
  end type matB

contains

  !----------------------------------------------------------------------!

  subroutine DQMC_B_Init(n, B, WS, Adj, t, mu, dtau)
    !
    ! Purpose
    ! =======
    !    This subroutine getnerates exponentional matrices of [-]dtau*T.
    !
    !         Ek  = exp(dtau*T)
    !         Eki = exp(-dtau*T)
    !
    !    The steps are
    !    1. Compute the Spectral Decomposition of T = USU^T
    !    2. Exponent the diagonal matrix S and -S
    !    3. Assemble B = U*exp(S)*U^T, Bi=U*exp(-S)*U^T
    !
    !    The parameters for checkboard method are also initialized
    !    in this subroutine. For the checkboard method, see [1] for
    !    more detail.
    !
    !    [1] W. Hanke and Yu.V. Kopaev, "Stable Numerical Simulation 
    !        of Models of Interating Electrons in Condensed-Matter Physics."
    !        Chapter 4. Elsevier Science Pub. 
    !
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)       :: n            ! Number of sites
    type(MatB), intent(inout) :: B            ! MatB
    type(WSpace), intent(in), target  :: WS   ! shared working space
    type(CCS), intent(in)     :: adj          ! adjacent info
    real(wp), intent(in)      :: t(:)         ! model parameter
    real(wp), intent(in)      :: mu(n), dtau  

    ! ... local scalars    ...
    integer  :: i, j                          ! iterator
    real(wp), pointer :: K(:,:)
    integer, pointer  :: start(:), r(:), A(:)

    ! ... Executable ...

    B%n    = n

    ! test if allocated
    allocate(B%B(n,n))
    allocate(B%Bi(n,n))

    B%name = "Dense MatB"

    ! Compute the B matrix
    K   => WS%R1

    K = ZERO
    start => Adj%cstart
    r     => Adj%row
    A     => Adj%A
    
    do i = 1, n
       do j = start(i), start(i+1)-1
          K(r(j),i) = t(A(j))*dtau
       end do
       K(i,i) = mu(i)*dtau
    end do

    call DQMC_B_ExpInit(B, K, WS)

  end subroutine DQMC_B_Init

  !----------------------------------------------------------------------!

  subroutine DQMC_B_ExpInit(B, K, WS)
    !
    ! Purpose
    ! =======
    !    This subroutine getnerates exponentional matrices of [-]dtau*T.
    !
    !         Ek  = exp(dtau*T)
    !         Eki = exp(-dtau*T)
    !
    !    The steps are
    !    1. Compute the Spectral Decomposition of T = USU^T
    !    2. Exponent the diagonal matrix S and -S
    !    3. Assemble B = U*exp(S)*U^T, Bi=U*exp(-S)*U^T
    !
    !    The parameters for checkboard method are also initialized
    !    in this subroutine. For the checkboard method, see [1] for
    !    more detail.
    !
    !    [1] W. Hanke and Yu.V. Kopaev, "Stable Numerical Simulation 
    !        of Models of Interating Electrons in Condensed-Matter Physics."
    !        Chapter 4. Elsevier Science Pub. 
    !
    !
    ! Arguments
    ! =========
    !
    type(MatB), intent(inout) :: B              ! MatB
    real(wp), intent(inout)   :: K(:,:)         ! model parameter
    type(WSpace), intent(in), target  :: WS     ! shared working space

    ! ... local scalars    ...
    integer  :: i, info, n       ! iterator
    real(wp), pointer :: W5(:), W4(:)

    ! ... Executable ...

    W5  => WS%R5
    W4  => WS%R7
    n   = B%n

    ! Compute the Spectral Decomposition of W1
    call dsyev('V', 'L', n, K, n, W5, W4, WS%lw(LA_SYEV), info)

    !! Check error message
    if (info .ne. 0) then
       call DQMC_Error("Error: from dsyev, dqmc_getek", info)
    end if
  
    ! Exponent the diagonal elements and scale the right eigenmatrix
    WS%R2 = K

    !! Scale the right eigenmatrix
    do i = 1, n
       call dscal(n, dexp(W5(i)), WS%R2(:,i), 1)
    end do

    ! Assemble B
    call dgemm('N', 'T', n, n, n, ONE, K, n, WS%R2, n, ZERO, B%B, n)


    ! For the inverse part, do the same thing\

    WS%R2 = K

    do i = 1, n
       call dscal(n, dexp(-W5(i)), WS%R2(:,i), 1)
    end do
    
    ! Assemble Bi
    call dgemm('N', 'T', n, n, n, ONE, K, n, WS%R2, n, ZERO, B%Bi, n)

  end subroutine DQMC_B_ExpInit

  !----------------------------------------------------------------------!


  subroutine DQMC_B_Free(B)
    !
    ! Purpose
    ! =======
    !    This subroutine frees memory of B.
    !
    ! Arguments
    ! =========
    !
    type(MatB), intent(inout)  :: B  ! MatB

    ! ... Executable ...
    
    deallocate(B%B, B%Bi)

  end subroutine DQMC_B_Free

  !----------------------------------------------------------------------!
  
  subroutine DQMC_MultB_Left(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix B_i from
    !    leftside.
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              !
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space   
    ! ... executable ...

    ! Multiply B from left-hand-side

    ! Copy M to C
    C = M

    ! M = B*C
    call dgemm('N', 'N', n, n, n, ONE, B%B, n, C, n, ZERO, M, n)

    ! M = Vi*M 
    call DQMC_ScaleRow(n, M, V_i)

  end subroutine DQMC_MultB_Left

  !-----------------------------------------------------------------------!

  subroutine DQMC_MultBi_Left(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix B_i from
    !    leftside.
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !    3. Flag lr is either 'l' or 'r'.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              ! MatB 
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space  
    ! ... executable ...

    ! Multiply B from left-hand-side

    ! M = inv(V_i)*M
    call DQMC_ScaleRowInv(n, M, V_i)

    ! Copy M to C
    C = M

    ! M = B*C
    call dgemm('N', 'N', n, n, n, ONE, B%Bi, n, C, n, ZERO, M, n)

  end subroutine DQMC_MultBi_Left

  !----------------------------------------------------------------------!
  
  subroutine DQMC_MultB_Right(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix B_i from
    !    rightside.
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !    3. Flag lr is either 'l' or 'r'.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              !
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space  
    ! ... executable ...

    ! Multiply B from right-hand-side
    ! M = M*V_i
    call DQMC_ScaleCol(n, M, V_i)
    
    ! Copy M to C
    C = M
    
    ! M = C*B
    call dgemm('N', 'N', n, n, n, ONE, C, n, B%B, n, ZERO, M, n)
    
  end subroutine DQMC_MultB_Right

  !-----------------------------------------------------------------------!

  subroutine DQMC_MultBi_Right(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix B_i from
    !    rightside.
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !    3. Flag lr is either 'l' or 'r'.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              ! MatB 
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space  
    ! ... executable ...

    ! Multiply B from right-hand-side
    ! Copy M to C
    C = M
    
    ! M = C*B
    call dgemm('N', 'N', n, n, n, ONE, C, n, B%Bi, n, ZERO, M, n)
    
    ! M = M*V_i
    call DQMC_ScaleColInv(n, M, V_i)
    
  end subroutine DQMC_MultBi_Right

  !-----------------------------------------------------------------------!

  subroutine DQMC_GetB(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine returns M = V_iB
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              !
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space  
    ! ... executable ...

    M = B%B
    ! M = Vi*M 
    call DQMC_ScaleRow(n, M, V_i)

  end subroutine DQMC_GetB

  !-----------------------------------------------------------------------!

  subroutine DQMC_GetBi(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine returns M = inv(B)inv(V_i)
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              ! MatB 
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         ! working space  
    ! ... executable ...

    M = B%Bi
    ! M = Vi*M 
    call DQMC_ScaleColInv(n, M, V_i)

  end subroutine DQMC_GetBi

  !-----------------------------------------------------------------------!

end module DQMC_MATB
