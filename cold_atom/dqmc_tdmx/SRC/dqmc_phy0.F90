module DQMC_Phy0
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_WSPACE
  use DQMC_STRUCT
  use LAPACK_MOD
  use BLAS_MOD

  implicit none 
  
  !
  ! This module contains data structure and subroutines for some 
  ! physics measurement on Hubbard's model, including
  ! 
  !     1.  Up spin occupancy
  !     2.  Down spin occupancy
  !     3.  Potential energy
  !     4.  Kinetic energy
  !     5.  Total engery
  !     6.  Density
  !     7.  XX ferromagnetic structure factor
  !     8.  ZZ ferromagnetic structure factor
  !     9.  XX antiferromagnetic structure factor
  !    10.  ZZ antiferromagnetic structure factor
  !    11.  RMS XX AF SF
  !    12.  RMS ZZ AF SF
  !    13.  Average sign
  !    14.  Average up sign
  !    15.  Average down sign
  !    16.  Equal time Green's function
  !    17.  Density-density correlation fn: (up-up)
  !    18.  Density-density correlation fn: (up-dn)
  !    19.  XX Spin correlation function
  !    20.  ZZ Spin correlation function
  !
  ! Measurement 1-15 are scalars. 16-20 are lists of length
  ! nClass, which are the number of distinct autocorrelation terms.
  ! (see DQMC_STRUCT for more details about nClass.)
  !
  ! List of subroutines
  ! ===================
  !    1. DQMC_Phy0_Init(P0, nClass, nBin, nHist)
  !    2. DQMC_Phy0_Normalize(P0, u)
  !    3. DQMC_Phy0_Print(P0, S, OPT)
  !    4. DQMC_Phy0_GetErr(P0)
  !    5. DQMC_Phy0_Dump(P0, idx, opt)
  !    6. DQMC_PHY0_Hist(n, nBin, over, under, H, list, GetIndex)
  !    7. DQMC_Meas0(n, P0, G_up, G_dn, mu, t, sgnup, sgndn, S, up, dn)
  !
  !    *** Most subroutines are for internal use only. User program 
  !        should not need any of them.
  !
  ! Data Structure and Parameters
  ! =============================
  !    The data type Phy0 is consisted of two parts: the measurements
  !    and the histograms. 
  !    
  !    The measurements are put into bins and then analysized in the end.
  !    To make data manipulation easy, all scalar variables are put into
  !    an array S and identified by the indeces, which are defined by the
  !    parameters below. Each column of S represents a bin.
  !    List variables, like Green's function, are put into separated 
  !    arrays, in which each column is also for one bin. 
  !        
  !    There are two special bins, averge and error, which are used to
  !    store the final averaged result and error. Their indeces are 
  !    specified by 'avg' and 'err' respectively.
  !
  !    The histogram part consists of three histograms, each having
  !    nHist+2 bins. They are histogram for up occupancy (Nup), 
  !    down occupancy (Ndn), and Nup*Ndn. The additional two bins
  !    are for data exceeds the range, whose indeces are specified
  !    by 'over' and 'under' respectively.
  !

  type Phy0
     ! Measurement part
     integer  :: nClass                     ! Number of distinct 
                                            ! autocorrelction terms
     integer  :: nBin                       ! Number of terms
     integer  :: nMeas                      ! Number of measurements
     integer  :: avg, err                   ! Index for average and error
                                            ! bins.
     integer  :: cnt                        ! Number of measurement for 
                                            ! current bin
     integer  :: idx                        ! current bin index
    
     ! Scalar array
     real(wp), pointer :: meas(:,:)        ! Scalar varaibles
     real(wp), pointer :: sign(:,:)        ! Scalar varaibles


     ! Array    
     real(wp), pointer :: G_fun(:,:)    ! Green's function
     real(wp), pointer :: SpinXX(:,:)   ! XX Spin correlation function
     real(wp), pointer :: SpinZZ(:,:)   ! ZZ Spin correlation function
     real(wp), pointer :: Den0(:,:)     ! Density-density correlation 
     real(wp), pointer :: Den1(:,:)     ! up-up (0) and up-dn (1) 

     ! working space
     real(wp), pointer :: up(:), dn(:)
     logical :: compSAF
     logical :: init

  end type Phy0

  ! Parameter for the index of scalar variables
  integer, parameter :: P0_NUP       = 1
  integer, parameter :: P0_NDN       = 2
  integer, parameter :: P0_NUD       = 3
  integer, parameter :: P0_KE        = 4
  integer, parameter :: P0_ENERGY    = 5
  integer, parameter :: P0_DENSITY   = 6

  integer, parameter :: P0_SFERRO    = 7
  integer, parameter :: P0_SFER2     = 8
  integer, parameter :: P0_SAF       = 9
  integer, parameter :: P0_SAFSQ     = 10
  integer, parameter :: P0_SAF2      = 11
  integer, parameter :: P0_SAF2SQ    = 12

  integer, parameter :: P0_N_NO_SAF  = 8
  integer, parameter :: P0_N         = 12

  integer, parameter :: P0_SGN       = 1
  integer, parameter :: P0_SGNUP     = 2
  integer, parameter :: P0_SGNDN     = 3


  ! Name of scalar variables
  character(len=*), parameter :: P0_STR(P0_N) = (/&
       "          Up spin occupancy : ", &
       "        Down spin occupancy : ", &
       "             <U*N_up*N_dn>  : ", &
       "             Kinetic energy : ", &
       "               Total Energy : ", &
       "                    Density : ", &
       "  XX Ferro structure factor : ", &
       "  ZZ Ferro structure factor : ", &
       "     XX AF structure factor : ", &
       "  Root Mean Square of XX AF : ", &
       "     ZZ AF structure factor : ", &
       "  Root Mean Square of ZZ AF : "/)

  character(len=*), parameter :: P0_SIGN_STR(3) = (/&
       "                   Avg sign : ", &
       "                Avg up sign : ", &
       "                Avg dn sign : "/)

contains

  ! Subroutines
  ! ==================================================================
  
  subroutine DQMC_Phy0_Init(P0, S, nBin, WS)
    !
    ! Purpose
    ! =======
    !    This subroutine initializes Phy0.
    !
    !  Pre-assumption
    ! ==============
    !    nClass, nBin and nHist are positive integers.
    !
    ! Arguments
    ! =========
    !
    type(Phy0), intent(inout) :: P0      ! Phy0 to be initialized
    type(Struct), intent(in)  :: S
    integer, intent(in)       :: nBin    ! No of bins
    type(WSpace), intent(in), target :: WS

    ! ... Local vars ...
    integer :: nClass

    ! ... Executable ...

    nClass     = S%nClass
    P0%nClass  = S%nClass
    P0%nBin    = nBin

    P0%avg     = nBin + 1
    P0%err     = nBin + 2
    P0%cnt     = 0
    P0%idx     = 1

    P0%compSAF = S%checklist(STRUCT_PHASE)
    ! count total number of data
    if (P0%compSAF) then
       P0%nMeas = P0_N
    else
       p0%nMeas = P0_N_NO_SAF
    end if

    ! Allocate storages
    !!   Allocate two additional bins for storing 
    !!   average value and error
    allocate(P0%meas(P0%nMeas, nBin+2))
    allocate(P0%sign(3, nBin+2))
    ! Initialize 
    P0%meas    = ZERO
    P0%sign    = ZERO

    allocate(P0%G_fun (nClass, nBin+2))
    allocate(P0%SpinXX(nClass, nBin+2))
    allocate(P0%SpinZZ(nClass, nBin+2))
    allocate(P0%Den0  (nClass, nBin+2))
    allocate(P0%Den1  (nClass, nBin+2))

    P0%G_fun   = ZERO
    P0%SpinXX  = ZERO
    P0%SpinZZ  = ZERO
    P0%Den0    = ZERO
    P0%Den1    = ZERO
    
    P0%up => WS%R5
    P0%dn => WS%R6

    P0%init = .true.
   end subroutine DQMC_Phy0_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy0_Free(P0)
    !
    ! Purpose
    ! =======
    !    This subroutine frees Phy0.
    !
    ! Arguments
    ! =========
    !
    type(Phy0), intent(inout) :: P0      ! Phy0 to be initialized

    ! ... Executable ...
    if (P0%init) then
       deallocate(P0%meas, P0%sign, P0%G_fun, P0%SpinXX, P0%SpinZZ)
       deallocate(P0%Den0, P0%Den1)
    end if

  end subroutine DQMC_Phy0_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy0_Avg(P0)
    !
    ! Purpose
    ! =======
    !    This subroutine averges the data with a bin.
    !    It also computes some measurements from others, such as
    !    density and energy.
    !
    !  Pre-assumption
    ! ==============
    !    idx is in the range of (1,nClass).
    !
    ! Arguments
    ! =========
    !
    type(Phy0), intent(inout) :: P0     ! Phy0
    
    ! ... local scalar ...
    real(wp) :: factor
    integer  :: idx, n

    ! ... Executable ...
    idx = P0%idx

    ! compute the normalization factor = 1/cnt
    if (P0%cnt .eq. 0) then
       call DQMC_Error("Phy0 normalize: cnt = 0", 0)
    end if
    factor = ONE/P0%cnt

    ! total occupancy = nup + ndn
    P0%meas(P0_DENSITY, idx) = P0%meas(P0_NUP, idx)+&
         P0%meas(P0_NDN, idx)
    ! total energy = kenitic energy + potential energy
    P0%meas(P0_ENERGY, idx) = P0%meas(P0_KE, idx)+&
         P0%meas(P0_NUD, idx)
    ! average
    P0%meas(:, idx) = P0%meas(:, idx) * factor
    P0%sign(:, idx) = P0%sign(:, idx) * factor
    
    if (P0%compSAF) then
       ! The sqaure terms
       P0%meas(P0_SAFSQ, idx)=sqrt(abs(P0%meas(P0_SAFSQ, idx)))
       P0%meas(P0_SAF2SQ,idx)=sqrt(abs(P0%meas(P0_SAF2SQ,idx)))
    end if

    ! This list terms
    n = P0%nClass
    call dscal(n, factor, P0%G_fun (1:n,idx), 1)
    call dscal(n, factor, P0%SpinXX(1:n,idx), 1)
    call dscal(n, factor, P0%SpinZZ(1:n,idx), 1)
    call dscal(n, factor, P0%Den0  (1:n,idx), 1)
    call dscal(n, factor, P0%Den1  (1:n,idx), 1)

    ! Change bin
    P0%idx = P0%idx + 1

    ! reset the counter
    p0%cnt = 0

  end subroutine DQMC_Phy0_Avg

  !--------------------------------------------------------------------!
  
  subroutine DQMC_Phy0_Print(P0, S, OPT)
    !
    ! Purpose
    ! =======
    !    This subroutine prints out the average and errors
    !    of measurements. Structure S will give labels for
    !    each autocorrelation terms.
    !
    !  Pre-assumption
    ! ===============
    !    OPT is a file handle
    !    DQMC_Phy0_GetErr was called.
    !
    ! Arguments
    ! =========
    !
    type(Phy0), intent(in)    :: P0   ! Phy0
    type(Struct), intent(in)  :: S    ! Underline lattice structure
    integer, intent(in)       :: OPT  ! Output file handle

    ! ... Local scalar ...
    integer :: nClass, avg, err

    ! ... Executable ...

    nClass = P0%nClass
    avg    = P0%avg
    err    = P0%err

    ! Scalar terms
    call DQMC_Print_RealArray(0, 3, "Sign of equal time measurements:", &
         P0_SIGN_STR, P0%sign(:,avg:avg), P0%sign(:,err:err), OPT)
    
    call DQMC_Print_RealArray(0, P0%nmeas, "Equal time measurements:", &
         P0_STR, P0%meas(:,avg:avg), P0%meas(:,err:err), OPT)

    ! Function terms
    call DQMC_Print_RealArray(0, nClass, "Equal time Green's function:", &
         S%clabel, P0%G_fun(:, avg:avg), P0%G_fun(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, &
         "Density-density correlation fn: (up-up)", &
         S%clabel, P0%Den0(:, avg:avg), P0%Den0(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, &
         "Density-density correlation fn: (up-dn)", &
         S%clabel, P0%Den1(:, avg:avg), P0%Den1(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, "XX Spin correlation function:", &
         S%clabel, P0%SpinXX(:, avg:avg), P0%SpinXX(:, err:err), OPT)
    
    call DQMC_Print_RealArray(0, nClass, "ZZ Spin correlation function:", &
         S%clabel, P0%SpinZZ(:, avg:avg), P0%SpinZZ(:, err:err), OPT)
    
  end subroutine DQMC_Phy0_Print

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy0_GetErr(P0)
    !
    ! Purpose
    ! =======
    !    This subroutine computes the average and errors
    !    of measurements. 
    !
    ! Argument
    ! ========
    !
    type(Phy0), intent(inout) :: P0

    ! ... Local Scalar ...
    integer  :: i, n
    integer  :: avg, err

    ! ... Local Array
    real(wp) :: sum_sgn, sgn(P0%nBin), y(P0%nBin), data(P0%nBin)
    
    ! ... Executable ...
    n   = P0%nBin
    avg = P0%avg
    err = P0%err

    ! Average sign, sign up, sign down
    do i = P0_SGNUP, P0_SGNDN
       data = P0%sign(i, 1:n)
       call DQMC_JackKnife(n, P0%sign(i, avg), P0%sign(i, err), data, &
            y, sgn, sum_sgn)
    end do
    
    data = P0%sign(P0_SGN, 1:n)
    call DQMC_JackKnife(n, P0%sign(P0_SGN, avg), P0%sign(P0_SGN, err), data, &
         y, sgn, sum_sgn)
    
    ! Average other single terms
    do i = 1, P0%nMeas
       data = P0%meas(i, 1:n)
       call DQMC_SignJackKnife(n, P0%meas(i, avg), P0%meas(i, err), data, &
            y, sgn, sum_sgn)
    end do
    
    ! Average Green function
    do i = 1, P0%nClass
       data =  P0%G_fun(i, 1:n)
       call DQMC_SignJackKnife(n, P0%G_fun(i, avg), P0%G_fun(i, err), &
            data, y, sgn, sum_sgn)
    end do
    
    ! Average correlated Density 
    do i = 1, P0%nClass
       data = P0%Den0(i, 1:n)
       call DQMC_SignJackKnife(n, P0%Den0(i, avg), P0%Den0(i, err), &
            data, y, sgn, sum_sgn)
    end do
    
    do i = 1, P0%nClass
       data = P0%Den1(i, 1:n)
       call DQMC_SignJackKnife(n, P0%Den1(i, avg), P0%Den1(i, err), &
            data, y, sgn, sum_sgn)
    end do
    
    ! Average spin
    do i = 1, P0%nClass
       data = P0%SpinXX(i, 1:n)
       call DQMC_SignJackKnife(n, P0%SpinXX(i, avg), P0%SpinXX(i, err), &
            data, y, sgn, sum_sgn)
    end do
    
    do i = 1, P0%nClass
       data = P0%SpinZZ(i, 1:n)
       call DQMC_SignJackKnife(n, P0%SpinZZ(i, avg), P0%SpinZZ(i, err), &
            data, y, sgn, sum_sgn)
    end do

  end subroutine DQMC_Phy0_GetErr

  !--------------------------------------------------------------------!

  subroutine DQMC_Phy0_Meas(n, P0, G_up, G_dn, U, mu, t, sgnup, sgndn, S)
    !
    ! Purpose
    ! =======
    !    This subroutine performs some physics measurement on
    !    Hubbard model.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)          :: n            ! Number of sites
    type(Phy0), intent(inout)    :: P0           ! Phy0
    real(wp), intent(in)         :: G_up(n,n)    ! Green's function
    real(wp), intent(in)         :: G_dn(n,n)    ! for spin up and down
    real(wp), intent(in)         :: sgnup, sgndn ! Sgn for det(G_up) det(G_dn)
    real(wp), intent(in)         :: mu(n), t(:)  ! Chemical and Kinetic para
    real(wp), intent(in)         :: U(:)         ! Chemical and Kinetic para
    type(Struct), intent(in)     :: S            ! Underline structure
    target :: S

    ! ... local scalar ...

    integer  :: i, j, k                          ! Loop iterator
    integer  :: tmp, idx, m                      ! Helper variables
    real(wp) :: sgn                        
    real(wp) :: var1, var2, var3          
    integer, pointer  :: start(:), r(:), A(:)

    ! ... executable ...

    idx = P0%idx
    tmp = P0%avg

    ! initialization
    ! Here we use avg bin as a temp variable 
    P0%meas(:,tmp)   = ZERO

    P0%G_fun(:,tmp)  = ZERO
    P0%Den0(:,tmp)   = ZERO
    P0%Den1(:,tmp)   = ZERO
    P0%SpinXX(:,tmp) = ZERO
    P0%SpinZZ(:,tmp) = ZERO
 

    ! Compute the site density for spin up and spin down
    do i = 1,n
       !======================================================!
       ! The density of electrons of spin up(dn) on site i    !
       ! is 1-G_up(i,i) (1-G_dn(i,i)).                        !
       ! nup (ndn) is the sum of all spin up (down) electrons.!
       !======================================================!
       P0%up(i)  = ONE - G_up(i,i)
       P0%dn(i)  = ONE - G_dn(i,i)
       P0%meas(P0_NUD, tmp) = P0%meas(P0_NUD, tmp)+ &
            P0%up(i)*P0%dn(i)*U(S%Map(i))
    end do

    P0%meas(P0_NUP, tmp) = sum(P0%up)
    P0%meas(P0_NDN, tmp) = sum(P0%dn)
    
    !=========================================!
    ! Kinetic energy = tt*sum_{ij}(G_ij+G_ji) !
    ! where site i and site j are neighbors   !
    !=========================================!
    ! set alias
    start => S%T%cstart
    r     => S%T%row
    A     => S%T%A
    
    ! loop all adj sites
    do i = 1, n  ! for each column
       do j = start(i), start(i+1)-1 ! for each nonzero elements
          P0%meas(P0_KE, tmp) =  P0%meas(P0_KE, tmp) + &
               t(A(j))*(G_up(r(j),i)+G_dn(r(j),i))
       end do
       P0%meas(P0_KE, tmp)  = P0%meas(P0_KE, tmp) - &
            mu(S%Map(i)) * (P0%up(i) + P0%dn(i))
    end do

    !=====================!
    ! Autocorelation term.!
    !=====================!
    if (P0%compSAF) then
       P0%meas(P0_SAF, tmp)  = TWO*n-P0%meas(P0_NUP, tmp)-&
            P0%meas(P0_NDN, tmp)
       P0%meas(P0_SAF2,tmp)  = P0%meas(P0_SAF, tmp)
    end if

    
    do i = 1,n
       do j = 1,n
          var1 = G_up(i,j) * G_up(j,i) + G_dn(i,j) * G_dn(j,i)
          var2 = -TWO * G_up(i,j) * G_dn(j,i)
          var3 = G_up(i,i) * G_up(j,j) + G_dn(i,i) * G_dn(j,j) - &
                 TWO * G_up(i,i) * G_dn(j,j) - var1
          
          ! k is the index
          k = S%D(i,j)
          P0%G_fun(k, tmp) = P0%G_fun(k, tmp) + G_up(i,j) + G_dn(i,j)
          P0%Den0(k, tmp)  = P0%Den0(k, tmp) + &
               P0%up(i)*P0%up(j) + P0%dn(i)*P0%dn(j) - var1
          P0%Den1(k, tmp)  = P0%Den1(k, tmp) + P0%up(i)*P0%dn(j)
          P0%SpinXX(k, tmp) = P0%SpinXX(k, tmp) + var2
          P0%SpinZZ(k, tmp) = P0%SpinZZ(k, tmp) + var3
          
          if (P0%compSAF) then
             var1 = S%P(i)*S%P(j)
             P0%meas(P0_SAF, tmp) = P0%meas(P0_SAF, tmp) + var1 * var2
             P0%meas(P0_SAF2,tmp) = P0%meas(P0_SAF2,tmp) + var1 * var3
          end if
       end do
       ! special case for (i,i)
       k = S%D(i,i)
       var1 =  G_up(i,i) + G_dn(i,i)
       P0%Den0(k, tmp)   = P0%Den0(k, tmp)   + var1
       P0%SpinXX(k, tmp) = P0%SpinXX(k, tmp) + var1
       P0%SpinZZ(k, tmp) = P0%SpinZZ(k, tmp) + var1
    end do
    
    P0%meas(P0_SFERRO, tmp) = sum(P0%SpinXX(:,tmp))
    P0%meas(P0_SFER2,  tmp) = sum(P0%SpinZZ(:,tmp))
    
    ! Average
    P0%meas(:,tmp) = P0%meas(:,tmp) / n
    do i = 1, P0%nClass
       P0%G_fun (i, tmp) = P0%G_fun (i, tmp) / S%F(i) * HALF
       P0%SpinXX(i, tmp) = P0%SpinXX(i, tmp) / S%F(i)
       P0%SpinZZ(i, tmp) = P0%SpinZZ(i, tmp) / S%F(i)
       P0%Den0  (i, tmp) = P0%Den0  (i, tmp) / S%F(i) * HALF
       P0%Den1  (i, tmp) = P0%Den1  (i, tmp) / S%F(i)
    end do
    
    if (P0%compSAF) then
       P0%meas(P0_SAFSQ, tmp) = P0%meas(P0_SAF, tmp) * P0%meas(P0_SAF, tmp)
       P0%meas(P0_SAF2SQ,tmp) = P0%meas(P0_SAF2,tmp) * P0%meas(P0_SAF2,tmp)
    end if

    ! Accumulate result to P0(:, idx)
    sgn = sgnup * sgndn
    P0%meas(:, idx) =  P0%meas(:, idx) + P0%meas(:, tmp) * sgn

    m = P0%nClass
    call daxpy(m, sgn, P0%G_fun (1:m,tmp), 1, P0%G_fun (1:m,idx), 1)
    call daxpy(m, sgn, P0%SpinXX(1:m,tmp), 1, P0%SpinXX(1:m,idx), 1)
    call daxpy(m, sgn, P0%SpinZZ(1:m,tmp), 1, P0%SpinZZ(1:m,idx), 1)
    call daxpy(m, sgn, P0%Den0  (1:m,tmp), 1, P0%Den0  (1:m,idx), 1)
    call daxpy(m, sgn, P0%Den1  (1:m,tmp), 1, P0%Den1  (1:m,idx), 1)

    P0%sign(P0_SGN,   idx) =  P0%sign(P0_SGN,   idx) + sgn
    P0%sign(P0_SGNUP, idx) =  P0%sign(P0_SGNUP, idx) + sgnup
    P0%sign(P0_SGNDN, idx) =  P0%sign(P0_SGNDN, idx) + sgndn
    P0%cnt = P0%cnt + 1
    
  end subroutine DQMC_Phy0_Meas

  !--------------------------------------------------------------------!

end module DQMC_Phy0
