module DQMC_TDM1
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_STRUCT
  use LAPACK_MOD
  use BLAS_MOD

  implicit none 
  
  !
  ! This module is designed for the computation of time dependent 
  ! measurement (TDM), which requires the unequal time Green's 
  ! function Gtau (up and down).
  !
  ! Two measurements are taken in this module: Green's function (G) and Chi
  ! (chi) function.  Both are stored in a 3 dimensional array. 
  ! The first dimension is for space, the second is for time, and the 
  ! last dimension is for bins.
  !
  type TDM1
     integer  :: n
     integer  :: nClass
     integer  :: L
     integer  :: tmp

     real(wp), pointer :: gnl(:,:,:)    ! Green's function, orignal
     real(wp), pointer :: chinl(:,:,:)  ! Chi function, orignal

     integer,  pointer :: D(:,:)
     integer,  pointer :: F(:)
  end type TDM1
  
contains

 !--------------------------------------------------------------------!
  
  subroutine DQMC_TDM1_Init(n, L, T1, nBin, S, tmp_idx)
    !
    ! Purpose
    ! =======
    !    This subroutine initializes TDM1. 
    !    Note the CS matrix is not initilized here.
    !    Need to initialize it separately.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1      ! time dependent measurement
    integer, intent(in)       :: n       ! No. of sites
    integer, intent(in)       :: L       ! No of time slice
    integer, intent(in)       :: nBin    ! No of Bins
    type(Struct), intent(in)  :: S
    target :: S
    integer, intent(in)       :: tmp_idx

    ! ... local variables ...
    integer     :: nClass

    ! ... Executable ...

    T1%n     = n
    nClass   = S%nClass
    T1%nClass= nClass
    T1%L     = L
    T1%D     => S%D
    T1%F     => S%F
    T1%tmp   = tmp_idx
    
    ! Allocate storages
    allocate(T1%gnl(nClass, L+1, nBin+2))
    allocate(T1%chinl(nClass, L+1, nBin+2))

    ! initialize values
    T1%gnl   = ZERO
    T1%chinl = ZERO
    
  end subroutine DQMC_TDM1_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Free(T1)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM1.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1      ! TDM to be freed

    ! ... Executable ...

    deallocate(T1%gnl, T1%chinl)

  end subroutine DQMC_TDM1_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Meas(T1, upt0, up0t, dnt0, dn0t, ti)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM1.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout)    :: T1
    real(wp), intent(in)         :: up0t(:,:), upt0(:,:)
    real(wp), intent(in)         :: dnt0(:,:), dn0t(:,:)
    integer, intent(in)          :: ti

    ! ... Local scalar ...

    integer  :: i, j, n, k
    integer,  pointer :: D(:,:)
    real(wp), pointer :: gnl(:), chinl(:)

    ! ... Executable ...

    ! Initialization
    n     =  T1%n
    D     => T1%D
    gnl   => T1%gnl(:,ti, T1%tmp)
    chinl => T1%chinl(:,ti, T1%tmp)

    ! Compute Green's function and Chi function
    do i = 1, n
       do j = 1, n
          ! k is the distance index of site i and site j
          k = D(i,j)
          gnl(k)  = gnl(k) + upt0(i,j) + dnt0(i,j)
          chinl(k)= chinl(k)-up0t(j,i)*dnt0(i,j) &
               - dn0t(j,i)*upt0(i,j)
       end do
    end do

  end subroutine DQMC_TDM1_Meas

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Acc(T1, sgn, idx)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM1.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout)    :: T1
    real(wp), intent(in)         :: sgn
    integer, intent(in)          :: idx

    ! ... Local scalar ...

    integer  :: i, j, L, nClass
    real(wp) :: factor
    real(wp), pointer :: gnl(:,:), chinl(:,:)
    integer, pointer  :: F(:)

    ! ... Executable ...

    ! Initialization
    L     = T1%L
    nClass= T1%nClass

    gnl   => T1%gnl(:,:,T1%tmp)
    chinl => T1%chinl(:,:,T1%tmp)
    F     => T1%F

    ! Average
    do j = 1, L
       do i = 1, nClass
          factor = sgn/F(i)
          T1%gnl(i,j,idx)   = T1%gnl(i,j,idx) + factor*gnl(i,j)*HALF
          T1%chinl(i,j,idx) = T1%chinl(i,j,idx) + factor*chinl(i,j)
       end do
    end do

    ! Special case for L+1
    T1%gnl(1,L+1,idx) = T1%gnl(1,L+1,idx) + &
                            (ONE - gnl(1,1)*HALF/F(1))*sgn
    T1%chinl(1,L+1,idx)=T1%chinl(1,L+1,idx)+sgn*chinl(1,1)/F(1)
    do i = 2, nClass
       factor = sgn/F(i)
       T1%gnl(i,L+1,idx)  = T1%gnl(i,L+1,idx) - factor*gnl(i,1)*HALF
       T1%chinl(i,L+1,idx)= T1%chinl(i,L+1,idx) + factor*chinl(i,1)
    end do
    
    ! Clean up
    T1%gnl(:,:,T1%tmp)   = ZERO
    T1%chinl(:,:,T1%tmp) = ZERO

  end subroutine DQMC_TDM1_Acc

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Avg(T1, idx, factor)
    !
    ! Purpose
    ! =======
    !    This subroutine averges the pair measurements.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1                 ! T1
    real(wp), intent(in)      :: factor
    integer, intent(in)       :: idx

    ! ... BLAS function ...
    real(wp), external :: ddot

    ! ... local scalar ...
    integer  :: nl

    ! ... Executable ...
    nl     =  T1%nClass*(T1%L+1)

    ! For Green's function and chi 

    ! Compute average on Green's function
    call dscal(nl, factor, T1%gnl(:,:,idx), 1)

    ! Compute average on Chi 
    call dscal(nl, factor, T1%chinl(:,:,idx), 1)

  end subroutine DQMC_TDM1_Avg

  !--------------------------------------------------------------------!

end module DQMC_TDM1
