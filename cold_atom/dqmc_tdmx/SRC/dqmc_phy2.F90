module DQMC_Phy2
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_STRUCT
  use DQMC_WSPACE
  use LAPACK_MOD
  use BLAS_MOD

  implicit none 
  
  !
  ! This module contains data structure and subroutines for 
  ! pair measurement.
  !
  ! List of subroutines
  ! ===================
  !    1. DQMC_Phy2_Init(P2, nBin, nb)
  !    2. DQMC_Phy2_Avg(P2, W, T, ldt)
  !    3. DQMC_Phy2_Print(P2, OPT)
  !    4. DQMC_Phy2_GetErr(P2)
  !    5. DQMC_Phy2_Meas(n, P2, G_up, G_dn, sgn, S)
  !
  ! Data Structure and Parameters
  ! =============================
  !    The data type Phy2 performs pair measurement. For each paired sites, 
  !    it computes the convolution of their correlation and their neighbors' 
  !    correlation. (CHECK) The results are accumulated into bins and then
  !    transformed by the 'wave' function. 
  !
  type Phy2
     ! Measurement of pair
     integer  :: nb                          ! Number of bonds
     integer  :: nWave                       ! Number of wave functions
     integer  :: nBin                        ! Number of bins
     integer  :: idx                         ! current bin index
     integer  :: cnt                         ! Number of measurement for 
                                             ! current bin 
     integer  :: avg, err                    ! Index for average and error
                                             ! bins.
     real(wp), pointer :: sgn(:)             ! sign  
     real(wp), pointer :: M1(:,:)            ! pair measurement
     real(wp), pointer :: M2(:,:)            ! accumulated pair measurement
     real(wp), pointer :: M3(:,:)            ! averaged pair measurement

     ! working space
     real(wp), pointer :: T(:,:)
     integer  :: ldt
     integer  :: nData
     logical  :: compute
  end type Phy2

contains

  ! Subroutines
  ! ==================================================================
  
  subroutine DQMC_Phy2_Init(P2, nBin, S, WS)
    !
    ! Purpose
    ! =======
    !    This subroutine initializes Phy2.
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(inout) :: P2          ! Phy2 to be initialized
    integer, intent(in)       :: nBin        ! No of bin
    type(Struct), intent(in)  :: S
    type(Wspace), intent(in), target :: WS   ! working space

    ! ... Local Vars ...
    integer :: nb, nWave

    ! ... Executable ...

    P2%compute = S%checklist(STRUCT_BOND)

    if (P2%compute) then
       nWave    = S%nWave
       nb     = S%n_b

       P2%nBin  = nBin
       P2%nb    = nb
       p2%nWave = nWave
       
       ! Allocate storages
       !! Allocate two additional bins for storing 
       !! average value and error
       allocate(P2%M1(nb,nb))
       allocate(P2%M2(nb,nb))
       allocate(P2%M3(nWave,nBin+2))
       allocate(P2%sgn(nBin+2))
       
       ! Initialize 
       P2%M2      = ZERO
       P2%sgn     = ZERO
       P2%avg     = nBin+1
       P2%err     = nBin+2
       P2%cnt     = 0
       P2%idx     = 1
       
       P2%T  =>  WS%R4
       P2%ldt     = size(WS%R4, 1)
       P2%nData   = nWave + 1
    else
       P2%nData = 0
    end if

  end subroutine DQMC_Phy2_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Free(P2)
    !
    ! Purpose
    ! =======
    !    This subroutine frees Phy2.
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(inout) :: P2      ! Phy0 to be initialized

    ! ... Executable ...

    if (P2%compute) then
       deallocate(P2%sgn, P2%M1, P2%M2, P2%M3)
    end if
 
  end subroutine DQMC_Phy2_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Meas(n, M1, M2, P2, Bond, G_up, G_dn, sgn)
    !
    ! Purpose
    ! =======
    !    This subroutine performs physics measurement on
    !    the neighbors of pairs of sites 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)          :: n             ! Number of sites
    real(wp), intent(inout)      :: M1(:, :)      ! temp variables
    real(wp), intent(inout)      :: M2(:, :)      !
    type(Phy2), intent(inout)    :: P2            ! Phy2
    type(CCS), intent(in)        :: Bond
    real(wp), intent(in)         :: G_up(n,n)     ! Green's function
    real(wp), intent(in)         :: G_dn(n,n)     ! for spin up and down
    real(wp), intent(in)         :: sgn           ! Sgnup*sgndn
 
    ! ... local scalar ...
    integer  :: nb
    real(wp) :: factor

    if (P2%compute) then

       ! Compute the pair measurement       
       call DQMC_Phy2_Pair(n, M1, G_up, G_dn, Bond)

       ! Averaging and accumulation
       factor = sgn/n
       nb = P2%nb
       call daxpy(nb*nb, factor, M1, 1, M2, 1)
       P2%sgn(P2%idx) = P2%sgn(P2%idx) + sgn
       
       ! increase counter
       p2%cnt = P2%cnt + 1
    end if

  end subroutine DQMC_Phy2_Meas

  !--------------------------------------------------------------------!
  
  subroutine DQMC_Phy2_Pair(n, M, up, dn, Bond)
    !
    ! Purpose
    ! =======
    !    This subroutine computes pair measurements.
    !    For each pair (i,j), compute the effect of 
    !    all their neighbor pairs.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)          :: n             ! Number of sites
    real(wp), intent(inout)      :: M(:, :)       ! pair measurement
                                                  ! (P_dd')
    real(wp), intent(in)         :: up(:,:)       ! Green's function
    real(wp), intent(in)         :: dn(:,:)       ! for spin up and down
    type(CCS), intent(in)        :: Bond          ! bond info
 
    ! ... local variables ...
    integer  :: i, j, inbr, jnbr, ni, nj, nci, ncj
    integer, pointer  :: start(:), r(:), A(:)

    ! ... Executable ...

    ! alias
    start => Bond%cstart
    r     => Bond%row
    A     => Bond%A
    
    ! clean up
    M = ZERO

    do i = 1, n
       do inbr = start(i), start(i+1)-1
          ni  = r(inbr)              ! bond of site i
          nci = A(inbr)              ! bond class of (i,ri)
          do j = 1, n
             do jnbr = start(j), start(j+1)-1
                nj  = r(jnbr)        ! bond of site j
                ncj = A(jnbr)        ! bond class of (j,rj)
                M(ncj, nci) = M(ncj, nci) + up(nj, ni)*dn(j, i)
             end do
          end do
       end do
    end do

  end subroutine DQMC_Phy2_Pair

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Avg(P2, W)
    !
    ! Purpose
    ! =======
    !    This subroutine averges the pair measurements.
    !    which is stored in P2%M2.
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(inout) :: P2                 ! phy2
    real(wp), intent(in)      :: W(P2%nb,P2%nwave)! wave

    ! ... local scalar ...
    real(wp) :: factor

    ! ... Executable ...

    if (P2%compute) then
       factor = ONE/P2%cnt       
       call DQMC_Wave_Avg(P2%nb, P2%nWave, W, P2%T, P2%M2, &
            P2%M3(1:P2%nWave,P2%idx), factor, P2%ldt)
       P2%sgn(P2%idx) = P2%sgn(P2%idx) * factor
    
       ! Reset counter and change bins
       P2%idx = P2%idx + 1
       P2%cnt = 0
       P2%M2  = ZERO
    end if

  end subroutine DQMC_Phy2_Avg

  !--------------------------------------------------------------------!

  subroutine DQMC_Wave_Avg(nb, nWave, W, T, I, O, factor, ldt)
    !
    ! Purpose
    ! =======
    !    This subroutine computes the waves from input data.
    !
    !    The averaging process runs as follows. 
    !    1. Let T = I*W
    !    2. Compute O(i) = W(:,i)'*T(:,i)
    !    3. Averaging O  = O*factor
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)     :: nB              ! number of bonds
    integer, intent(in)     :: nWave           ! number of waves
    real(wp), intent(in)    :: W(nb, nWave)    ! wave matrix
    real(wp), intent(in)    :: T(ldt, nWave)   ! temp matrix
    real(wp), intent(in)    :: I(nb, nb)       ! input matrix
    real(wp), intent(inout) :: O(:)            ! output matrix
    real(wp), intent(in)    :: factor
    integer, intent(in)     :: ldt

    ! ... Local Variables ...
    integer :: j

    ! ... BLAS function ...
    real(wp), external :: ddot

    ! ... Executable ....

    call dgemm('N','N', nb, nWave, nb, factor, I, nb, W, nb, &
         ZERO, T, ldt)
    do j = 1, nWave
       O(j) = ddot(nb, W(1:nb,j), 1, T(1:nb,j), 1)
    end do

  end subroutine DQMC_Wave_Avg


  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_GetErr(P2)
    !
    ! Purpose
    ! =======
    !    This subroutine computes the average and errors
    !    of measurements. 
    !
    ! Argument
    ! ========
    type(Phy2), intent(inout) :: P2

    ! ... Local Scalar ...
    integer  :: i
    integer  :: n, avg, err

    ! ... Local Array
    real(wp) :: sgn(P2%nBin) , sum_sgn, y(P2%nBin), data(P2%nBin)
    
    ! ... Executable ...

    if (P2%compute) then
       n   = P2%nBin
       avg = P2%avg
       err = P2%err
       
       ! Averge the other terms
       data = P2%sgn(1:n)
       call DQMC_JackKnife(n, P2%sgn(avg), P2%sgn(err), data, &
            y, sgn, sum_sgn)
       
       do i = 1, P2%nWave
          data = P2%M3(i, 1:n)
          call DQMC_SignJackKnife(n, P2%M3(i, avg), P2%M3(i, err), &
               data, y, sgn, sum_sgn)
       end do
       
    end if
    
  end subroutine DQMC_Phy2_GetErr

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Print(P2, wlabel, OPT)
    !
    ! Purpose
    ! =======
    !    This subroutine prints out the average and errors
    !    of pair measurements. The first three terms only.
    !
    !  Pre-assumption
    ! ===============
    !    OPT is a file handle
    !    DQMC_Phy2_GetErr was called.
    !    nWave >= 3
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(in)    :: P2   ! Phy2
    character(*), intent(in)  :: wlabel(:)
    integer, intent(in)       :: OPT  ! Output file handle

    ! ... Local var
    integer  :: avg, err
    
    ! ... Executable ...

    if (P2%compute) then
       avg = P2%avg
       err = P2%err
       call DQMC_Print_RealArray(0, P2%nWave, "Pair measurement:", &
            wlabel, P2%M3(:,avg:avg), P2%M3(:,err:err), OPT)
    end if
    
  end subroutine DQMC_Phy2_Print

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Pack(P2, darray, start)
    !
    ! Purpose
    ! =======
    !    Packs all the data of P0 into a data array.
    !    darray is an (k, nbin) array, where k is the total number
    !    of measurement of entire simulation.
    !    start is the starting index for P0.
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(in)    :: P2           ! Phy0
    real, intent(inout)       :: darray(:,:)
    integer, intent(in)       :: start        ! starting point of darray
  
    ! ... local scalar ...
    integer :: i1, i2, nBin, nWave

    ! ... executable ...

    if (P2%compute) then
       nWave = P2%nWave
       nBin = P2%nBin
       darray(start, 1:nBin) = real(P2%sgn(1:nBin))
       
       i1   = start+1
       i2   = start + nWave
       darray(i1:i2, 1:nBin) = real(P2%M3(1:nWave, 1:nBin))
       
    end if

  end subroutine DQMC_Phy2_Pack

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy2_Unpack(P2, darray, start)
    !
    ! Purpose
    ! =======
    !    Unpacks data from darray to P0.
    !    darray is an (k, nbin) array, where k is the total number
    !    of measurement of entire simulation.
    !    start is the starting index for P0.
    !
    ! Arguments
    ! =========
    !
    type(Phy2), intent(inout) :: P2           ! Phy0
    real, intent(in)          :: darray(:,:)
    integer, intent(in)       :: start        ! starting point of darray
  
    ! ... local scalar ...
    integer :: i1, i2, nBin, nWave

    ! ... executable ...

    if (P2%compute) then
       nWave = P2%nWave
       nBin = P2%nBin
       P2%sgn(1:nBin) = darray(start, 1:nBin)
       
       i1   = start+1
       i2   = start + nWave
       P2%M3(1:nWave, 1:nBin) = darray(i1:i2, 1:nBin)
    end if

  end subroutine DQMC_Phy2_Unpack

  !--------------------------------------------------------------------!

end module DQMC_Phy2
