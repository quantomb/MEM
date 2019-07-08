module DQMC_CKB
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_WSPACE
  use BLAS_MOD
  use LAPACK_MOD

  implicit none 

  ! 
  ! This module defines the type and subroutines for propagator B.
  ! B is defined by the checkerboard method. See [1]
  !
  ! For a 2-dimensional model
  !
  ! B = exp(K) = exp(K_{XEVEN})exp(K_{YEVEN})exp(K_{YODD})exp(K_{XODD})
  !            = B_4B_3B_2B_1
  !            
  !    [1] W. Hanke and Yu.V. Kopaev, "Stable Numerical Simulation 
  !        of Models of Interating Electrons in Condensed-Matter Physics."
  !        Chapter 4. Elsevier Science Pub. 
  !
  !  Data Type
  !  =========
  !
  type matB
     integer  :: n                      ! dim of B
     integer  :: m                      ! number of neighbors of lattice
     integer, pointer :: A(:,:)         ! Adjacency info
     real(wp),pointer :: factor(:)      ! temp variable in multiplication
     real(wp) :: soc, f1, f2            ! parameters for checkerboard method
     character(12) :: name
  end type matB

contains


  !-------------------------------------------------------------------------!

  subroutine DQMC_B_Init(n, B, WS, Adj, BS, nBand, t, mu, dtau)
    !
    ! Purpose
    ! =======
    !    This subroutine initiliazes the data type of Green function.
    !
    ! Pre-assumption
    ! ==============
    !    This module only used for one band Hubbard's model.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)       :: n                   ! Number of sites
    type(MatB), intent(inout) :: B                   ! MatB
    type(WSpace), intent(in)  :: WS                  ! shared working space
    integer, intent(in),target       :: adj(:,:)     ! adjacent info
    integer, intent(in)       :: nBand               ! Number of band
    integer, intent(in)       :: BS(:,:)             ! band size
    real(wp), intent(in)      :: t(:)                ! model parameter
    real(wp), intent(in)      :: mu(:), dtau  

    ! ... local scalars    ...
    real(wp) :: tdtau, ch
    
    ! ... Executable ...

    B%n    = n
    B%m    = BS(1,1)
    B%name = "Checkerboard"

    tdtau   = t(1)*dtau
    ch      = cosh(tdtau)
    B%soc   = sinh(tdtau)/ch
    if (mu(1) .eq. ZERO) then
       B%f1 = ch**B%m
       B%f2 = B%f1
    else
       B%f1 = exp(dtau*mu(1))*(ch**B%m)
       B%f2 = (ch**B%m)/exp(dtau*mu(1))
    end if
    B%A     => Adj
    
    call DQMC_Reshape(n, B%factor)

  end subroutine DQMC_B_Init

  !-------------------------------------------------------------------------!

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

    nullify(B%A)
    call DQMC_Free(B%factor)

  end subroutine DQMC_B_Free


  !----------------------------------------------------------------------!
  
  subroutine DQMC_MultB_Left(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply B from
    !    leftside.
    !
    !      B_i = V_i*B
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
    integer, intent(in)          :: n              ! The order of matrix M  
    real(wp), intent(inout)      :: M(n,n)         !  
    type(MatB), intent(in)       :: B              !
    real(wp), intent(in)         :: V_i(n)         !
    real(wp), intent(inout)      :: C(n,n)         ! working space     

    ! ... Local variables ...
    integer  :: i, j, k, k2, m1
    real(wp) :: soc, f1
    integer,  pointer :: A(:,:) 
    real(wp), pointer :: factor(:)

    ! ... executable ...

    A      => B%A
    soc    = B%soc
    factor => B%factor
    m1     = B%m
    f1     = B%f1

    ! Combin V_i with the coefficient of B
    do i = 1,n
       factor(i) = f1*V_i(i)
    end do
    
    ! Multiply B from left-hand-side

    do j = 1, n ! do it column by column         
       ! Preassumption: B%m must be even
       do k = 1, m1, 2

          ! To avoid matrix copy, we use the temporary matrix C
          ! store the immediate result of M.
          ! C = B_{i}*M, and then 
          ! M = B_{i+1}*C = B_{i+1}*B_{i}*M
          do i = 1, n
             C(i,1) = M(i,j) + soc*M(A(i,k),j)
          end do
          
          k2 = k + 1
          do i = 1, n
             M(i,j) = C(i,1) + soc*C(A(i,k2),1)
          end do
          
       end do
       
       ! Multiply V          
       do i = 1, n
          M(i,j) = factor(i) * M(i,j)
       end do
    end do

  end subroutine DQMC_MultB_Left

  !----------------------------------------------------------------------!
  
  subroutine DQMC_MultB_Right(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply B from
    !    righthand side.
    !
    !      B_i = V_i*B
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
    integer, intent(in)          :: n              ! The order of matrix M  
    real(wp), intent(inout)      :: M(n,n)         !  
    type(MatB), intent(in)       :: B              !
    real(wp), intent(in)         :: V_i(n)         !
    real(wp), intent(inout)      :: C(n,n)         !    

    ! ... Local variables ...
    integer  :: i, j, k, k2
    real(wp) :: soc
    integer,  pointer :: A(:,:) 
    real(wp), pointer ::factor(:)

    ! ... executable ...

    A      => B%A
    soc    = B%soc
    factor => B%factor
    
    ! Combin V_i with the coefficient of B
    do i = 1,n
       factor(i) = B%f1*V_i(i)
    end do

    ! Multiply B from right-hand-side

    ! Multiply V
    call DQMC_ScaleCol(n, M, factor)

    ! To optimize the cache performance.
    ! we transpose M first, and then transpose it back.later
    call DQMC_Trans(n, C, M)
    
    ! M = M*B = (B'*M')' = B_m'*...*B_2'*B_1'*M
    do j = 1, n
       ! Preassumption: B%m must be even
       do k = B%m, 1, -2
          
          ! To avoid matrix copy, we use the temporary matrix C
          ! store the immediate result of M.
          ! C = B_{i}*M, and then 
          ! M = B_{i+1}*C = B_{i+1}*B_{i}*M
          do i = 1, n
             M(i,1) = C(i,j) + soc*C(A(i,k),j)
          end do
          
          k2 = k - 1
          do i = 1, n
             C(i,j) = M(i,1) + soc*M(A(i,k2),1)
          end do

       end do

    end do
    ! transpose back
    call DQMC_Trans(n, M, C)

  end subroutine DQMC_MultB_Right

  !--------------------------------------------------------------------------!
  
  subroutine DQMC_MultBi_Left(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply inv(B) from
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
    integer, intent(in)          :: n              ! The order of matrix M  
    real(wp), intent(inout)      :: M(n,n)         !  
    type(MatB), intent(in)       :: B              !
    real(wp), intent(in)         :: V_i(n)         !
    real(wp), intent(inout)      :: C(n,n)         !
    
    ! ... Local variables ...
    integer  :: i, j, k, k2
    real(wp) :: soc 
    integer,  pointer :: A(:,:) 
    real(wp), pointer :: factor(:)

    ! ... executable ...

    A      => B%A
    soc    = B%soc
    factor => B%factor

    ! Combin V_i with the coefficient of B
    do i = 1, n
       factor(i) = B%f2/V_i(i)
    end do

    ! Multiply B from left-hand-side
    ! M = inv(B_1)*...inv(B_m)*M
    do j = 1, n
       
       ! Multiply V
       do i = 1, n
          M(i,j) = factor(i) * M(i,j)
       end do

       ! Preassumption: B%m must be even
       do k = B%m, 1, -2
          
          do i = 1, n
             C(i,1) = M(i,j) - soc*M(A(i,k),j)
          end do
          
          k2 = k - 1
          do i = 1, n
             M(i,j) = C(i,1) - soc*C(A(i,k2),1)
          end do
          
       end do
    end do
    
  end subroutine DQMC_MultBi_Left

  !--------------------------------------------------------------------------!
  
  subroutine DQMC_MultBi_Right(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply inv(B) from
    !    right hand side.
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)          :: n              ! The order of matrix M  
    real(wp), intent(inout)      :: M(n,n)         !  
    type(MatB), intent(in)       :: B              !
    real(wp), intent(in)         :: V_i(n)         !
    real(wp), intent(inout)      :: C(n,n)         !
    
    ! ... Local variables ...
    integer  :: i, j, k, k2
    real(wp) :: soc 
    integer,  pointer :: A(:,:) 
    real(wp), pointer :: factor(:)

    ! ... executable ...

    A      => B%A
    soc    = B%soc
    factor => B%factor

    ! Combin V_i with the coefficient of B
    do i = 1, n
       factor(i) = B%f2/V_i(i)
    end do

    ! To optimize the cache performance.
    ! we transpose M first, and then transpose it back later       
    call DQMC_Trans(n, C, M)
    
    ! M = M*inv(B) = (inv(B)'*M')' = (inv(B_m)*...*inv(B_1)*M')'
    do j = 1, n
       ! Preassumption: B%m must be even
       do k = 1, B%m, 2
          
          ! To avoid matrix copy, we use the temporary matrix C
          ! store the immediate result of M.
          ! C = B_{i}*M, and then 
          ! M = B_{i+1}*C = B_{i+1}*B_{i}*M
          do i = 1, n
             M(i,1) = C(i,j) - soc*C(A(i,k),j)
          end do
          
          k2 = k + 1
          do i = 1, n
             C(i,j) = M(i,1) - soc*M(A(i,k2),1)
          end do
       end do

    end do
    
    ! transpose back
    call DQMC_Trans(n, M, C)

    ! Multiply V
    call DQMC_ScaleCol(n, M, factor)
    
  end subroutine DQMC_MultBi_Right

  !-------------------------------------------------------------------------!

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
    real(wp), intent(inout)  :: C(n,n)         !
    ! ... Executable ...

    call DQMC_Eye(n, M)
    call DQMC_MultB_Left(n, M, B, V_i, C)

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
    real(wp), intent(inout)  :: C(n,n)         !
    ! ... Executable ...

    call DQMC_Eye(n, M)
    call DQMC_MultBi_Left(n, M, B, V_i, C)

  end subroutine DQMC_GetBi

  !-----------------------------------------------------------------------!
  
end module DQMC_CKB
