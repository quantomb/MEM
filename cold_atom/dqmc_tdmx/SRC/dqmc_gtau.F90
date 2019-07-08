module DQMC_GTAU
#include "dqmc_include.h"

  use DQMC_UTIL
  use LAPACK_MOD
  use BLAS_MOD
  use DQMC_WSPACE
  use _DQMC_MATB
  use DQMC_SEQB
  use DQMC_GFUN

  implicit none 
  
  !
  ! This module is designed for the computation of time dependent 
  ! measurement (TDM), which requires the unequal time Green's 
  ! function Gup_tau (Gdn_tau).
  ! Mathematically, Gup_tau(Gdn_tau) is the inverse of the following
  ! matrix.
  !
  !                  [  I                   B_L] 
  !                  [-B_1   I                 ]
  !              M = [     -B_2  I             ]
  !                  [          ...   ...      ]
  !                  [             -B_{L-1}  I ]
  !                   
  !
  ! The current implementation only considers to return the ith row
  ! of Gup_tau (Gdn_tau). Let Gup_ij be the (i,j) block matrix of 
  ! Gup_tau.
  !
  !     Gup_ii = inv(I+B_iB_{i-1}...B_1B_L...B_{i+1})   
  !     Gup_ij = -Gup_ii(B_iB_{i-1}...B_1B_L...B_{j+1})  for j = i...L  
  !     Gup_ij =  Gup_ii(B_iB_{i-1}...B_1B_L...B_{j+1})  for j = 1...i-1  
  !
  ! The corresponding Gdn_ij has the same structure.
  !
  ! [1] Z. Bai, W.Chen, R. Scalettar, I. Yamazaki, "Lecture Notes 
  !     on Advances of Numerical Methods for Hubbard Quantum Monte
  !     Carlo Simulation." 
  !
  type Gtau
     integer  :: n                           ! dimension of G_j
     integer  :: L                           ! number of block columns
     integer  :: ii, ib                      ! index of the block
     integer  :: sfc                         ! safe count
     integer  :: nWrap                       ! safe wrapping
     integer  :: which                       ! which part of gtau
                                             ! should be perfomed

     ! working space for constructing unequal time Green's function
     real(wp), pointer :: upt0(:,:)          ! used for UDT decomp
     real(wp), pointer :: up0t(:,:)
     real(wp), pointer :: dnt0(:,:)
     real(wp), pointer :: dn0t(:,:)
     real(wp), pointer :: U(:,:)             ! used in UDT decomp
     real(wp), pointer :: T(:,:)
     real(wp), pointer :: D(:) 
     real(wp), pointer :: v1(:), v2(:)       ! used in G computing
     real(wp), pointer :: v3(:), v4(:)
     real(wp), pointer :: W1(:,:)

     type(SeqB) :: SB1, SB2
     type(MatB), pointer :: B
     
     real(wp) :: sgnup, sgndn

  end type Gtau
  
  integer, parameter :: TAU_T0   = 0   ! Column
  integer, parameter :: TAU_BOTH = 1   ! column and row
  integer, parameter :: TAU_0T   = 2   ! ROW

contains

  ! Subroutines
  ! ==================================================================
  
  subroutine DQMC_Gtau_Init(n, L, which, nOrth, nWrap, tau, B, WS)
    !
    ! Purpose
    ! =======
    !    This subroutine initializes Phy2.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)       :: n, L
    integer, intent(in)       :: which   ! 
    integer, intent(in)       :: nOrth
    integer, intent(in)       :: nWrap
    type(Gtau), intent(inout) :: tau     ! time dependent measurement
    type(WSpace), intent(in)  :: WS
    type(MatB), target, intent(in)    :: B

    ! ... Executable ...

    tau%n     = n
    tau%L     = L
    tau%which = which
    tau%nWrap = nWrap
    tau%sfc   = 0
    tau%ii    = 0
    tau%ib    = 0

    ! Allocate storages
    allocate(tau%upt0(n,n))
    allocate(tau%dnt0(n,n))
    allocate(tau%up0t(n,n))
    allocate(tau%dn0t(n,n))
    allocate(tau%U(n,n))
    allocate(tau%T(n,n))
    allocate(tau%D(n))
    allocate(tau%v1(n))
    allocate(tau%v2(n))
    allocate(tau%v3(n))
    allocate(tau%v4(n))

    call DQMC_SeqB_Init(n, L, nOrth, B, tau%SB1, WS)
    call DQMC_SeqB_Init2(n, L, nOrth, B, &
         tau%SB2, tau%U, tau%D, tau%T, WS)

    tau%SB2%piv => WS%I2
    tau%W1 => WS%R1
    tau%B  => B

  end subroutine DQMC_Gtau_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_Gtau_Free(tau)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM.
    !
    ! Arguments
    ! =========
    !
    type(Gtau), intent(inout) :: tau      ! TDM to be initialized

    ! ... Executable ...

    deallocate(tau%upt0, tau%up0t, tau%dnt0, tau%dn0t)
    deallocate(tau%U, tau%D, tau%T, tau%v1, tau%v2, tau%v3, tau%v4)

    call DQMC_SeqB_Free(tau%SB1)
    call DQMC_SeqB_Free(tau%SB2)
 
  end subroutine DQMC_Gtau_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_MakeGtau(tau, G_up, G_dn, ii, ib)
    !
    ! Purpose
    ! =======
    !    This subroutine generates Gtau.
    !
    ! Arguments
    ! =========
    !
    type(Gtau), intent(inout)    :: tau
    type(G_fun), intent(inout)   :: G_up, G_dn
    integer, intent(in)          :: ii, ib     ! the index of block matrix
    
    ! ... local scalar
    integer  :: n, idx, i
    real(wp) :: det
    logical  :: recompute

    ! ... executable ...
    !
    !  meaning of indices
    !     ii: the ii-th block row or column
    !     ib: the block offset
    !
    !     id: 

    ! initialization
    n = tau%n
    recompute = .false.
    
    ! start computation
    if (ib .gt. 0) then ! Off diagoanl block

       ! If the current Gtau is not a neighbor of requested one,
       ! set count=0 for recomputing.
       if (ii .eq. tau%ii .and. ib .eq. tau%ib+1) then
          ! recude safe count
          tau%sfc = tau%sfc - 1
       else
          ! recompute if cannot use update
          tau%sfc = 0
       end if

       ! idx is the absolute row/column index
       idx = ii+ib
       if (idx .gt. tau%L) then
          idx = idx - tau%L
       end if


       ! compute Gtau
       if (tau%sfc .ne. 0) then 
          ! Within a toerlable range, the Gtau can be updated
          ! by simple matrix multiplication.
          ! see [1] for more details.
          
          ! use update
          if (tau%which .le. TAU_BOTH) then
             
             ! column             
             call DQMC_MultB_Left(n, tau%upt0, tau%B, G_up%V(:,idx), &
                  tau%W1)
             call DQMC_MultB_Left(n, tau%dnt0, tau%B, G_dn%V(:,idx), &
                  tau%W1)
          end if
          
          if (tau%which .ge. TAU_BOTH) then
             ! row
             call DQMC_MultBi_Right(n, tau%up0t, tau%B, G_up%V(:,idx), &
                  tau%W1)
             call DQMC_MultBi_Right(n, tau%dn0t, tau%B, G_dn%V(:,idx), &
                  tau%W1)
          end if
       else
          ! we need to change index, becuse the index system used in get
          ! gtau is absolute. (ib,jb)
          ! Recompute Gtau from scratch
          
          call DQMC_GetGtau(ii, idx, tau%upt0, tau%up0t, G_up%V, tau)
          call DQMC_GetGtau(ii, idx, tau%dnt0, tau%dn0t, G_dn%V, tau)
          recompute = .true.
          tau%sfc = tau%nWrap

       end if

    else ! The initial ii block
       ! Construct G from getG function
       ! In this case, upt0 = -up0t and dnt0 = -dn0t       
       call DQMC_ComputeG(ii, n, tau%sgnup, tau%upt0, G_up%V, tau%SB1, &
            G_up%pvt, .false., det)   
       call DQMC_ComputeG(ii, n, tau%sgndn, tau%dnt0, G_dn%V, tau%SB1, &
            G_up%pvt, .false., det) 
       
       ! if compute both, then use up0t as companion matrices
       if (tau%which .eq. TAU_BOTH) then
          tau%up0t = tau%upt0
          tau%dn0t = tau%dnt0
          do i = 1, n
             tau%up0t(i,i) = tau%up0t(i,i) - ONE
             tau%dn0t(i,i) = tau%dn0t(i,i) - ONE
          end do
       end if
       
       ! set the counter
       tau%sfc = tau%nWrap
       recompute = .true.
    end if
    
    ! determine the sign
    if ((recompute .and. ii+ib.gt.tau%L).or.(ii+ib.eq.tau%L+1))then
       tau%upt0 = -tau%upt0
       tau%dnt0 = -tau%dnt0
       tau%up0t = -tau%up0t
       tau%dn0t = -tau%dn0t
    end if
    
    tau%ii = ii
    tau%ib = ib

  end subroutine DQMC_MakeGtau

  !--------------------------------------------------------------------------!
  ! The following subroutines are used for time dependent Green's function.  !
  !--------------------------------------------------------------------------!
  subroutine DQMC_GetGtau(ib, jb, G_ij, G_ji, V, tau)
    !
    ! Purpose
    ! =======
    !    
    !    This subroutine computes the (i,j) submatrix of Gtau if which
    !    equals to'R'ow or 'B'oth. and computes the (j,i) submatrix of 
    !    Gtau if which equals to 'C'olumn or 'B'oth 
    !
    ! Mathematically, Gtau(Gdn_tau) is the inverse of 
    !
    !                  [  I                   B_1] 
    !                  [-B_2   I                 ]
    !              M = [     -B_3  I             ]
    !                  [          ...   ...      ]
    !                  [             -B_{L}    I ]
    !                   
    !
    ! The (i,j) submatrix of Gtau is given as
    !
    !     G_ii =    inv(I+B_iB_{i-1}...B_1B_L...B_{i+1})   
    !     G_ij = -Gtau_ii(B_iB_{i-1}...B_1B_L...B_{j+1})  for j = i+1...L  
    !     G_ij =  Gtau_ii(B_iB_{i-1}...B_{j+1})           for j = 1...i-1 
    !
    ! In general, we can write G_ij as
    !
    !         G_ij = (+/-) inv(I+A_1A_2)A_1
    !              = (+/-) inv(inv(A_1)+A_2)
    ! where
    !
    !          A_1 = B_{i}...B_{j+1} and
    !          A_2 = B_{j}...B_{i+1}
    !
    ! The following procedure compute G_ij in a stable way
    ! 
    ! 1. Perform UDT decomposition on inv(A_1) and A_2
    !    
    !       inv(A_1) = U_1D_1T_1
    !           A_2  = U_2D_2T_2
    !
    !    See the DQMC_UDTD in DQMC_B.F90 for detail of UDT decomposition.
    !
    ! 2. Decompose D_1 = barD_1*hatD_1
    !              D_2 = barD_2*hatD_2
    !    where
    !           barD_1(i,i) = max(1, D_1(i,i)) and
    !           hatD_1(i,i) = min(1, D_1(i,i))
    !
    ! 3. Compute
    !
    !    C = hatD_2*T_2*inv(T_1)*inv(barD_1)+inv(barD_2)*inv(U_2)*U_1*hatD_1
    !    
    ! 4. Assemble G as 
    !
    !    G = inv(T_1)*inv(barD_1)*inv(C)*inv(barD_2)*inv(U_2)
    !      = inv(T_1)*inv(D_2T_2inv(T_1)+inv(U_2)U_1D_1)*inv(U_2)
    !      = inv(U_1D_1T_1+U_2D_2T_2)
    !      = inv(inv(A_1)+A_2)
    !      = inv(I+A_1A_2)A_1
    !      = inv(I+B_{i}...B_1*B_l...B_{i-1})B_i...B_{j+1}
    !
    ! Matrix G_ji has very similar structure with G_ij.
    !
    !     G_jj =    inv(I+B_jB_{j-1}...B_1B_L...B_{j+1})   
    !     G_ji = -Gtau_jj(B_jB_{j-1}...B_1B_L...B_{i+1})  for i = j+1...L  
    !     G_ji =  Gtau_jj(B_jB_{j-1}...B_{i+1})           for i = 1...j-1 
    !
    ! For a fixed i and j,
    !
    !         G_ji = (+/-) inv(I+A_2A_1)A_2
    !              = (+/-) inv(inv(A_2)+A_1)
    ! 
    ! where A_1 and A_2 are as defined before.
    ! Therefore, 
    !
    !         G_ji = inv(inv(U_1D_1T_1)+inv(U_2D_2T_2))
    !              = inv(inv(T_1)inv(D_1)inv(U_1)+inv(T_2)inv(D_2)inv(U_2))
    !              = U_2*inv(inv(D_1)inv(U_1)U_2+T_1*inv(T_2)*inv(D_2))*T_1
    !
    ! The same trick of bar and hat is also applied to inv(D_1) and inv(D_2).
    !        
    !         G_ji = U_2*inv(barD_2)*inv(...)*inv(barD_1)*T_1
    !
    ! where (...) = hatD_1*inv(U_1)U_2*inv(barD_2)+
    !               inv(barD_1)T_1*inv(T_2)hatD_2
    !
    ! NOTE: the hatD_1, barD_1, hatD_2, and barD_2 here are different from
    !       the previous ones.
    !
    ! See working notes for more detail.
    !
    ! Arguments
    ! =========
    !
    integer,  intent(in)       :: ib, jb        ! block indices
    real(wp), intent(inout)    :: G_ij(:,:)     ! submatrix of Gtau
    real(wp), intent(inout)    :: G_ji(:,:)     ! submatrix of Gtau
    real(wp), intent(inout)    :: V(:,:)        ! HSF
    type(Gtau), intent(inout)  :: tau

    ! ... local scalars    ...
    integer :: info           ! parameters for lapack's sub
    integer :: i              ! iterator
    integer :: n

    real(wp), pointer :: U1(:,:)       ! 
    real(wp), pointer :: D1(:)         ! 
    real(wp), pointer :: T1(:,:)       ! 
    real(wp), pointer :: U2(:,:)       ! 
    real(wp), pointer :: D2(:)         ! 
    real(wp), pointer :: T2(:,:)       ! 
    real(wp), pointer :: W1(:,:)       ! working space
    real(wp), pointer :: W2(:,:)       !
    real(wp), pointer :: rw(:)         ! working space
    integer,  pointer :: lw(:)         !
    integer,  pointer :: pvt1(:)       !
    integer,  pointer :: pvt2(:)       !

    real(wp), pointer :: bar1i(:)      ! 
    real(wp), pointer :: bar2i(:)      ! 
    real(wp), pointer :: hat1(:)       ! 
    real(wp), pointer :: hat2(:)       ! 
    

    ! ... Executable ...

    ! STEP 0. Initialization
    n = tau%n
    bar1i => tau%v1
    bar2i => tau%v2
    hat1  => tau%v3
    hat2  => tau%v4
    
    U1 => tau%SB1%U
    D1 => tau%SB1%D
    T1 => tau%SB1%T

    U2 => tau%SB2%U
    D2 => tau%SB2%D
    T2 => tau%SB2%T

    W1 => tau%SB1%W1
    W2 => tau%SB1%W2
    rw => tau%SB1%rw
    lw => tau%SB1%lw
    pvt1 => tau%SB1%piv
    pvt2 => tau%SB2%piv

    info = 0

    ! STEP 1. Cmpute UDT decomposition of 
    !         inv(A_1) = inv(B_{i}...B_{j+1})
    !         and A_2  = B_j...B_{i+1}.
    ! ==========================================
    ! W1, W2, rw, lwork, tau, pvt1 can be reused.

    call DQMC_SeqMultB (ib, jb+1, tau%SB1, V)
    call DQMC_SeqMultBi(jb, ib+1, tau%SB2, V)
    
    if (tau%which .eq. TAU_T0 .or. tau%which .eq. TAU_BOTH) then
       !
       ! STEP 2.  D_1 = inv(barD_1)*hatD_1
       !          D_2 = inv(barD_2)*hatD_2
       ! ==================================
       
       do i = 1, n
          bar1i(i) = ONE / max(ONE, abs(D1(i)))
          hat1(i)  = D1(i) * bar1i(i)
          bar2i(i) = ONE / max(ONE, abs(D2(i)))
          hat2(i)  = D2(i) * bar2i(i)
       end do

       !   
       ! STEP 3. Compute C = hatD_2*T_2*inv(T_1)*inv(barD_1)+
       !                     inv(barD_2)*inv(U_2)*U_1*hatD_1
       ! =======================================================   
       
       !! Compute  T_2*inv(T_1)
       ! copy T_1 to W_2, because we may need T_1 later
       call dcopy(n*n,T1,1,W2,1)
       
       ! W_1 = T_2'
       call DQMC_trans(n, W1, T2)

       ! W_1 = inv(W_2')*W_1 = inv(T_1')*T_2'
       call dgetrf(n, n, W2, n, pvt1, info)
       call dgetrs('T', n, n, W2, n, pvt1, W1, n, info)
       if (info .ne. 0) then
          call DQMC_Error("Error: dgetrs(1) in dqmc_getgtau.", info)
       end if

       ! T_2 = transpose(W_1) = transpose(inv(T_1')*T_2') = T_2*inv(T_1)
       call DQMC_trans(n, T2, W1)
       
       if (tau%which.eq.TAU_T0) then
          ! if only Row is computed, then T1 is not reference, reuse it 
          call dcopy(n*n,G_ij,1,T1,1)
       end if
       
       ! U_1 = G_ij = U_2'*U_1
       ! ** G_ij here is used as a temp variable
       call dgemm('T', 'N', n, n, n, ONE, U2, n, U1, n, ZERO, G_ij, n)
       call dcopy(n*n,G_ij,1,U1,1)
       
       !! *** We need to keep T2 and U1 for later use.
       
       ! compute U_1 = barD_2*U_2'*U_1*hatD_1
       call DQMC_ScaleRow(n, G_ij, bar2i)
       call DQMC_ScaleCol(n, G_ij, hat1)
       
       ! compute W_1 = hatD_2*T_2*inv(T_1)*barD_1
       call dcopy(n*n,T2,1,W1,1)
       call DQMC_ScaleRow(n, W1, hat2)
       call DQMC_ScaleCol(n, W1, bar1i)
       
       ! W_1 = W_1 + G_ij
       call daxpy(n*n, ONE, G_ij, 1, W1, 1)
       
       !   
       ! STEP 4. Compute inv(T_1)*inv(barD_1)*inv(C)*inv(barD_2)*inv(U_2)
       ! =================================================================   
       
       ! Let G_ij = inv(barD_2) * inv(U2)
       call DQMC_trans(n, G_ij, U2)
       call DQMC_ScaleRow(n, G_ij, bar2i)

       ! G_ij = inv(W_1)*inv(barD_2)*inv(U_2)
       call dgesv(n, n, W1, n, pvt2, G_ij, n, info)

       if (info .ne. 0) then
          call DQMC_Error("Error: dgesv(2) in dqmc_getgtau.", info)
       end if
       
       ! G_ij = inv(barD_1)*G_ij = inv(barD_1)*inv(W_1)*inv(barD_2)*inv(U_2)
       call DQMC_ScaleRow(n, G_ij, bar1i)
       
       ! G_ij = inv(T_1)*G_ij
       !      = inv(T_1)*inv(barD_1)*inv(C)*inv(barD_2)*inv(U_2)
       call dgetrs('N', n, n, W2, n, pvt1, G_ij, n, info)
       if (info .ne. 0) then
          call DQMC_Error("Error: dgetrs(1) in dqmc_getgtau.", info)
       end if

    end if

    !
    ! Compute G_ji, repeat step 2, 3, 4 for Gji
    ! ==========================================

    if (tau%which.eq.TAU_0T .or. tau%which .eq. TAU_BOTH) then       
       !
       ! STEP 5.  inv(D_1) = barD_1*hatD_1
       !          inv(D_2) = barD_2*hatD_2
       ! ======================================
       !
       do i = 1, n
          if (D1(i) .eq. ZERO) then
             call DQMC_Error("Error: in dqmc_getgtau, D1(i)=0.0, i=", i)
          end if
          D1(i) = ONE / D1(i)
          bar1i(i) = ONE / max(ONE, D1(i))
          hat1(i) = D1(i) * bar1i(i)

          if (D2(i) .eq. ZERO) then
             call DQMC_Error("Error: in dqmc_getgtau, D2(i)=0.0, i=", i)
          end if
          D2(i) = ONE / D2(i)
          bar2i(i) = ONE / max(ONE, D2(i))
          hat2(i) = D2(i) * bar2i(i)
       end do
       
       !   
       ! STEP 6. Compute G_ji = hatD_1*inv(U_1)U_2*inv(barD_2)+
       !                        inv(barD_1)T_1*inv(T_2)hatD_2
       ! =======================================================   
       if (tau%which .eq. TAU_BOTH) then
          ! Previously, T_2 = T_2*inv(T_1)
          !             U_1 = inv(U_2)*U_1
          ! Therefore, we only need to invert them.
          
          ! first, compute inv(barD_1)T_1*inv(T_2)hatD_2          
          call dgetrf(n, n, T2, n, pvt1, info)
          if (info .ne. 0) then
             call DQMC_Error("Error: dgetrf(1) in dqmc_getgtau.", info)
          end if
          call dgetri(n, T2, n, pvt1, rw, lw(LA_GETRI), info)
          
          ! W1 = U1' = inv(inv(U_2)*U_1) = inv(U_1)*U_2
          call DQMC_Trans(n, W1, U1)
           
       else
          ! No previous computed results. Compute them from scratch.
          !
          ! (1) Compute T_1*inv(T_2) 
          !     Let W_1 = T_1'
          call DQMC_trans(n, W1, T1)

          !     W_1 = inv(T_2')*W_1 = inv(T_2')*T_1'
          call dgetrf(n, n, T2, n, pvt1, info)
          call dgetrs('T', n, n, T2, n, pvt1, W1, n, info)
          if (info .ne. 0) then
             call DQMC_Error("Error: dgetrs(1) in dqmc_getgtau.", info)
          end if
          !     T_2 = W_1' = (inv(T_2')*T_1')' = T_1*inv(T_2)
          call DQMC_trans(n, T2, W1)
          
          ! (2) Compute W_1 = U_2'*U_1
          call dgemm('T', 'N', n, n, n, ONE, U1, n, U2, n, ZERO, W1, n)
       end if


       ! Compute inv(barD_1)T_1*inv(T_2)hatD_2
       call DQMC_ScaleRow(n, T2, bar1i)
       call DQMC_ScaleCol(n, T2, hat2)
          
       ! Compute hatD_1*inv(U_1)U_2*inv(barD_2)
       call DQMC_ScaleRow(n, W1, hat1)
       call DQMC_ScaleCol(n, W1, bar2i)
       

       ! W1 = W1 + T2
       call daxpy(n*n, ONE, T2, 1, W1, 1)

       !
       ! STEP 7. Compute U_2*inv(barD_2)*inv(...)*inv(barD_1)*T_1
       ! =========================================================
       !
       call DQMC_ScaleRow(n, T1, bar1i)
       call dgesv(n, n, W1, n, pvt1, T1, n, info)
       if (info .ne. 0) then
          call DQMC_Error("Error: dgesv(3) in dqmc_getgtau.", info)
       end if
       
       ! inv(barD_2)*inv(...)*inv(barD_1)*T_1
       call DQMC_ScaleRow(n, T1, bar2i)

       ! copy the previous result
       call dcopy(n*n,G_ij,1,W2,1)

       ! multiply -U2, the sign is negative
       call dgemm('N', 'N', n, n, n, -ONE, U2, n, T1, n, ZERO, G_ji, n)

    end if
    
  end subroutine DQMC_GetGtau

  !--------------------------------------------------------------------!
  ! This is only for test purpose                                      !
  !--------------------------------------------------------------------!

  subroutine DQMC_Gtau_Big(tau, Aup, Adn, G_up, G_dn)
    !
    ! This solves the gtau explicitly by Lapack.
    !

    type(Gtau), intent(inout)    :: tau
    type(G_fun), intent(inout)   :: G_up, G_dn
    real(wp), intent(inout), target :: Aup(:,:), Adn(:,:)
    
    ! ... Local var ...
    integer, parameter :: BS = 64
    integer  :: i, n, L, nL, lwork, info, ipiv(tau%n*tau%L)
    real(wp) :: work(BS*tau%n*tau%L)
    real(wp), pointer :: up(:,:), dn(:,:)
    
    ! ... Executable ...
    
    n   = tau%n
    L   = tau%L
    nL  = n*L
    Aup = ZERO
    Adn = ZERO
    ipiv = 0
    info = 0
    work = ZERO

    ! making A
    do i = 1, nL
       Aup(i,i) = ONE
       Adn(i,i) = ONE
    end do

    up => Aup(1:n, nL-n+1:nL)
    dn => Adn(1:n, nL-n+1:nL)

    call DQMC_GetB(n, up, tau%B, G_up%V(:,1), tau%W1)
    call DQMC_GetB(n, dn, tau%B, G_dn%V(:,1), tau%W1)

    do i = 1, L-1
       up =>  Aup(i*n+1:(i+1)*n, (i-1)*n+1:i*n)
       dn =>  Adn(i*n+1:(i+1)*n, (i-1)*n+1:i*n)
       call DQMC_GetB(n, up, tau%B, G_up%V(:,i+1), tau%W1)
       call DQMC_GetB(n, dn, tau%B, G_dn%V(:,i+1), tau%W1)
       up = -up
       dn = -dn
    end do

    ! invert A
    lwork = BS*nL
    call DGETRF(nL, nL, Aup, nL, ipiv, info)
    call DGETRI(nL, Aup, nL, ipiv, work, lwork, info)
    call DGETRF(nL, nL, Adn, nL, ipiv, info)
    call DGETRI(nL, Adn, nL, ipiv, work, lwork, info)
  end subroutine DQMC_Gtau_Big

end module DQMC_GTAU
