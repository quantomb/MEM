module DQMC_TDMX_RT
#include "dqmc_include.h"

  use DQMC_util
  use DQMC_Cfg
  use DQMC_STRUCT
  use LAPACK_MOD
  use BLAS_MOD
  use DQMC_Gtau

  implicit none 
  
  !
  ! This module is designed for the computation of time dependent 
  ! measurement (TDM), which requires the unequal time Green's 
  ! function Gtau (up and down).
  !
  ! Two measurements are taken in this module: 
  ! dynmod_a(il,idx): <c^+_{i,s}(t)c_{neibor(i')s}(0)>.<c_{neibor(i),s}(0)c^+_{i',s}(t)>
  ! dynmod_b(il,idx): <c^+_{i,s}(t)c_{neibor(i')s}(t)>.<c^+_{i',s}(0)c_{neibor(i),s}(0)>
  ! The first dimension il is for time, the second is for bins
  !
  type TDMX_RT
  ! time dependent Green's function matrix
  ! tau_eqt: equal-time Green's function
  ! tau_uneqt: unequal_time Green's function
     type(gtau) :: tau_eqt,tau_uneqt
  ! n: number of total size
     integer  :: n
  ! nClass: number of unique distance for time-dependent quantity
     integer  :: nClass
  ! L: number of time slices
     integer  :: L
  ! itvl: fine grained time slice number
     integer  :: itvl
  ! dtau: interval of time slice
     real(wp) :: dtau
  ! alphax : ratio of delta U/ delta J
     real(wp) :: alphax
  ! avgsgn: average sign of each data bin
     real(wp) :: sgn
  ! cnt: number of measurement in each data bin
     integer  :: cnt
  ! temperary index and average index
     integer  :: tmp,avg
  ! NADJ: number of adjacencies per site
     integer :: NADJ
  ! neibro(i,j) : adjuacency table of site j
     integer,  pointer :: neibor(:,:)
  ! unique distance table
     integer,  pointer :: D(:,:)
     integer,  pointer :: F(:)
  ! Green's function
     real(wp), pointer :: gnl(:,:,:)
  ! dynmod_a is used to measure \chi(\tau) in dynamical modulation
  ! index of dynmod_tkk(dr,dT,(avg or tmp index))
     real(wp), pointer :: dynmod_tkk(:,:,:)
     real(wp), pointer :: dynmod_tkd(:,:,:)
     real(wp), pointer :: dynmod_tdd(:,:,:)
     real(wp), pointer :: dynmod_t(:,:,:)
  ! dynmod_w0 is the temperature-dependent susceptibility
  ! dynmod_w(type,index), type=KK, KD+DK, DD, EE
     real(wp), pointer :: dynmod_w0(:,:)
  ! dynmod_e(index): the value of operator 
  ! dynmod_e(1,:) --> <K>; dynmod_e(2,:) --> <D>
     real(wp), pointer :: dynmod_e(:,:)

  ! FTM is Fourier Transform matrix
  ! FTM(T_i,iw_n)
     real(wp), pointer :: FTM(:,:)
     real(wp), pointer :: CS(:,:)
     logical :: FTspace   
  end type TDMX_RT
  
contains

 !--------------------------------------------------------------------!
  
  subroutine DQMC_TDMX_Init(TX, S, WS, B, cfg, alphax)
    !
    ! Purpose
    ! =======
    !    This subroutine initializes TDMX. 
    !    Note the CS matrix is not initilized here.
    !    Need to initialize it separately.
    !
    ! Arguments
    ! =========
    !
    type(TDMX_RT), intent(inout) :: TX   ! time dependent measurement
    type(Struct), intent(in)  :: S       ! Structure of Lattice
    target :: S
    type(MatB), intent(in)    :: B       ! B matrix
    type(WSpace), intent(in)  :: WS
    type(config), intent(in)  :: cfg
    real(wp), intent(in) :: alphax

    ! ... local variables ...
    integer :: n,L,nOrth,nWrap,itvl,size_bin
    integer :: i,j,nClass, longL
    integer, pointer :: start(:),r(:),A(:)
    integer, pointer :: find(:)
    real(wp) :: dtau

    ! ... Executable ...
    call CFG_Get(cfg, "n",n)
    call CFG_Get(cfg, "L",L)
    call CFG_Get(cfg, "north",nOrth)
    call CFG_Get(cfg, "nwrap",nWrap)
    call CFG_Get(cfg, "dtau",dtau)
    call CFG_Get(cfg, "nitvl",itvl)

    TX%L     = L
    TX%n     = n
    nClass   = S%nClass
    TX%nClass= nClass
    TX%tmp   = 1
    TX%avg   = 2
    TX%dtau  = dtau
    TX%alphax= alphax
    TX%itvl  = itvl
    TX%cnt   = 0
    TX%sgn   = ZERO
    longL = (L-1)*itvl + 1

    ! find neighborhood matrix
    ! transfer adjacencies table
    start => S%T%cstart
    r     => S%T%row
    A     => S%T%A
    ! determine number of adjacencies
    TX%NADJ = S%T%nnz/S%T%n
    ! allocate neibor
    allocate(TX%neibor(TX%NADJ,n))
    allocate(find(n))
    ! initial array for number of adjacency to be found 
    do i = 1, n
       find(i) = 0
    enddo
    do i = 1, n
       do j = start(i), start(i+1)-1
          ! find one adjacency
          find(r(j))=find(r(j)) + 1
          TX%neibor(find(r(j)),r(j)) = i
       enddo
    enddo
    ! test neibor table
    !do i = 1, n
    !   do j = 1, TX%NADJ
    !      write(*,*) j, i, TX%neibor(j,i)
    !   enddo
    !enddo

    !Allocate unique distance table
    !allocate(TX%D(nClass,nClass))
    !allocate(TX%F(nClass))
    ! Allocate storages
    allocate(TX%gnl(nClass,L,2))
    allocate(TX%dynmod_tkk(nClass,L,2))
    allocate(TX%dynmod_tkd(nClass,L,2))
    allocate(TX%dynmod_tdd(nClass,L,2))
    allocate(TX%dynmod_t(nClass,L,2))
    allocate(TX%dynmod_w0(4,2))
    allocate(TX%dynmod_e(2,2))
    allocate(TX%FTM(longL,2*L))

    ! initialize values
    TX%gnl = ZERO
    TX%dynmod_tkk = ZERO
    TX%dynmod_tkd = ZERO
    TX%dynmod_tdd = ZERO
    TX%dynmod_t = ZERO
    TX%dynmod_w0 = ZERO
    TX%dynmod_e = ZERO
    TX%FTM      = ZERO
    TX%D => S%D 
    TX%F => S%F  
    ! Form array of Fourier Transformation
    call DQMC_Make_FTMX(L,itvl,dtau,TX,ONE)

    ! initial gtau!!!
    call DQMC_Gtau_Init(n,L,TAU_BOTH,nOrth,nWrap,TX%tau_eqt,  B,WS)
    call DQMC_Gtau_Init(n,L,TAU_BOTH,nOrth,nWrap,TX%tau_uneqt,B,WS)

    ! allocate for FT
    TX%FTspace = S%checklist(STRUCT_FT)
    if(TX%FTspace) then
       TX%CS => S%FT
    else
       nullify(TX%CS)
    endif

    ! deallocate find
    deallocate(find)
    
  end subroutine DQMC_TDMX_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_TDMX_Free(TX)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM1.
    !
    ! Arguments
    ! =========
    !
    type(TDMX_RT), intent(inout) :: TX      ! TDM to be freed

    ! ... Executable ...
    ! Free Gtau
    call DQMC_Gtau_Free(TX%tau_eqt)
    call DQMC_Gtau_Free(TX%tau_uneqt)

    deallocate(TX%gnl)
    deallocate(TX%dynmod_tkk, TX%dynmod_tkd, TX%dynmod_t)
    deallocate(TX%dynmod_tdd)
    deallocate(TX%dynmod_w0)
    deallocate(TX%dynmod_e)
    deallocate(TX%neibor, TX%FTM)


  end subroutine DQMC_TDMX_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_TDMX_Meas(tm, G_up, G_dn)
    !
    ! Purpose
    ! =======
    !
    ! Arguments
    ! =========
    !
    type(TDMX_RT), intent(inout)   :: tm
    type(G_fun), intent(inout)  :: G_up, G_dn


    ! ... Local var ...
    integer :: i,n, L
    real(wp):: sgn
    real(wp),pointer :: up0t_eq(:,:), upt0_eq(:,:), dn0t_eq(:,:), dnt0_eq(:,:)
    real(wp),pointer :: up0t_un(:,:), upt0_un(:,:), dn0t_un(:,:), dnt0_un(:,:)
    ! special case for t=0
    ! G_kk(1:n) : G_{k,k}(0,0) i.e. equal-time Green's function for site k
    ! G_kl(1,n) : Sum of Green's function for neighbour k=x=l, G_{k,l}(0,0)
    real(wp),allocatable :: Gt0up(:,:)
    real(wp),allocatable :: Gt0dn(:,:)
    real(wp) :: dtau
    integer  :: itvl

    ! ... executable ...
    L = tm%L
    n = tm%n
    allocate(Gt0up(n,n))
    allocate(Gt0dn(n,n))
    Gt0up = ZERO 
    Gt0dn = ZERO
    itvl = tm%itvl
    dtau = tm%dtau


    upt0_eq => tm%tau_eqt%upt0
    up0t_eq => tm%tau_eqt%up0t
    dnt0_eq => tm%tau_eqt%dnt0
    dn0t_eq => tm%tau_eqt%dn0t

    upt0_un => tm%tau_uneqt%upt0
    up0t_un => tm%tau_uneqt%up0t
    dnt0_un => tm%tau_uneqt%dnt0
    dn0t_un => tm%tau_uneqt%dn0t

    !!! Special case for T_i = 0
    call DQMC_MakeGtau(tm%tau_eqt,G_up,G_dn,1,0)
    call DQMC_TDMX_GNL(tm,upt0_eq,up0t_eq,dnt0_eq,dn0t_eq,1)
    call DQMC_TDMX_EQT(tm,upt0_eq,up0t_eq,dnt0_eq,dn0t_eq,Gt0up,Gt0dn)

    !!! for case of T_i =x= 0
    do i = 2, L

    ! make equal-time Green's function for <c^+(t_i)c(t_i)>, t_i=i
       call DQMC_MakeGtau(tm%tau_eqt,   G_up, G_dn, i, 0)

    ! make unequal-time Green's function for <c^+(t_i)c(0)>, t_i=i
       call DQMC_MakeGtau(tm%tau_uneqt, G_up, G_dn, 1, i-1)

    ! measure unequal-time Green's function for <c^+(t_i)c(0)>, t_i=i
       call DQMC_TDMX_GNL(tm,upt0_un,up0t_un,dnt0_un,dn0t_un,i)

    ! measure susceptibilities
       call DQMC_TDMX_MeasALL(tm, upt0_eq, up0t_eq, dnt0_eq, dn0t_eq,&
            upt0_un, up0t_un, dnt0_un, dn0t_un, Gt0up, Gt0dn, i)
    end do

    sgn = G_up%sgn*G_dn%sgn
    ! Normalize
    call DQMC_TDMX_FTM(dtau,tm)
    ! Accumulate X(T) and X(iw_n)
    call DQMC_TDMX_Acc(tm,sgn)

    tm%sgn = tm%sgn + sgn
    tm%cnt = tm%cnt + 1

    deallocate(Gt0up,Gt0dn)

  end subroutine DQMC_TDMX_Meas
  !--------------------------------------------------------------------!

  subroutine DQMC_TDMX_GNL(TX, upt0, up0t, dnt0, dn0t, ti)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM1.
    !
    ! Arguments
    ! =========
    !
    type(TDMX_RT), intent(inout)    :: TX
    real(wp), intent(in)         :: up0t(:,:), upt0(:,:)
    real(wp), intent(in)         :: dnt0(:,:), dn0t(:,:)
    integer, intent(in)          :: ti

    ! ... Local scalar ...

    integer  :: i, j, n, k
    integer,  pointer :: D(:,:)
    real(wp), pointer :: gnl(:) !, chinl(:)

    ! ... Executable ...

    ! Initialization
    n     =  TX%n
    D     => TX%D
    gnl   => TX%gnl(:,ti, TX%tmp)
    !chinl => TX%chinl(:,ti, T1%tmp)

    ! Compute Green's function and Chi function
    do i = 1, n
       do j = 1, n
          ! k is the distance index of site i and site j
          k = D(i,j)
          gnl(k)  = gnl(k) + upt0(i,j) + dnt0(i,j)
          !chinl(k)= chinl(k)-up0t(j,i)*dnt0(i,j) &
          !     - dn0t(j,i)*upt0(i,j)
       end do
    end do

  end subroutine DQMC_TDMX_GNL
  !--------------------------------------------------------------------!

  subroutine DQMC_TDMX_MeasALL(TX,upt0_eq,up0t_eq,dnt0_eq,dn0t_eq,&
             upt0_un,up0t_un,dnt0_un,dn0t_un,Gt0up,Gt0dn,ti)
    !
    ! Purpose
    ! =======
    ! This subroutine measures the susceptibility
    ! <K(T_i)K(0)> - alphx*<K(T_i)D(0)> - alphax*<D(T_i)K(0)> 
    ! + alphax^2*<D(T_i)D(0)> 
    ! at time slice T_i
    ! 
    ! Arguments
    ! =========
    ! upt0_eq, up0t_eq, dnt0_eq, dn0t_eq: equal time GF at T=T_i
    ! upt0_un, up0t_un, dnt0_un, dn0t_un: unequal time GF at T=T_i
    ! Gsum_kl(1): sum of spin-up G_kl(0,0) when k =x= l
    ! Gsum_kl(2): sum of spin-down G_kl(0,0) when k =x= l
    ! G_kk(1,i): spin-up G_kk(0,0) at site i
    ! G_kk(2,i): spin-down G_kk(0,0) at site i
    !
    type(TDMX_RT), intent(inout)    :: TX
    real(wp),   intent(in)       :: up0t_eq(:,:), upt0_eq(:,:)
    real(wp),   intent(in)       :: dnt0_eq(:,:), dn0t_eq(:,:)
    real(wp),   intent(in)       :: up0t_un(:,:), upt0_un(:,:)
    real(wp),   intent(in)       :: dnt0_un(:,:), dn0t_un(:,:)
    real(wp),   intent(in)       :: Gt0up(:,:),Gt0dn(:,:)
    integer,    intent(in)       :: ti

    ! ... Local scalar ...

    integer  :: i, j, n, k, NADJ, nClass
    integer  :: ib, jb, i_neibor, j_neibor
    integer,  pointer :: D(:,:)
    ! neighborhood matrix
    integer,  pointer :: neibor(:,:)
    !it is not a array since it is just the single-valued function of \tau
    real(wp), pointer :: dynmod_tkk(:), dynmod_tkd(:)
    real(wp), pointer :: dynmod_t(:), dynmod_tdd(:)
    real(wp), allocatable :: value_kk(:),value_kd(:)
    real(wp), allocatable :: value_dk(:),value_dd(:)
    real(wp) :: alphax,tmp1,tmp2

    ! ... Executable ...

    ! Initialization
    n        =  TX%n
    nClass   =  TX%nClass
    neibor   => TX%neibor
    D        => TX%D
    dynmod_tkk => TX%dynmod_tkk(:,ti,TX%tmp)
    dynmod_tkd => TX%dynmod_tkd(:,ti,TX%tmp)
    dynmod_tdd => TX%dynmod_tdd(:,ti,TX%tmp)
    dynmod_t   => TX%dynmod_t(:,ti,TX%tmp)
    NADJ     =  TX%NADJ
    alphax   =  TX%alphax

    allocate(value_kk(nClass))
    allocate(value_kd(nClass))
    allocate(value_dk(nClass))
    allocate(value_dd(nClass))

    value_kk = ZERO
    value_kd = ZERO
    value_dk = ZERO
    value_dd = ZERO

    !write(*,*) '-------------------------------------------------------'
    !write(*,*) 'In T_i=',ti,' case'
    !write(*,*) 'upt0_un(1,2)=',upt0_un(1,2),' dnt0_un(1,2)=',dnt0_un(1,2)
    !write(*,*) 'up0t_un(1,2)=',up0t_un(1,2),' dn0t_un(1,2)=',dn0t_un(1,2)
    !write(*,*) 'upt0_un(2,1)=',upt0_un(2,1),' dnt0_un(2,1)=',dnt0_un(2,1)
    !write(*,*) 'up0t_un(2,1)=',up0t_un(2,1),' dn0t_un(2,1)=',dn0t_un(2,1)
    !pause
    ! Compute dynamic susceptibility for dynamical modulation
    ! <K(T_i)K(0)>
    do i = 1, n
    do ib = 1, NADJ, 2
       i_neibor = neibor(ib,i)
       do j = 1, n
       do jb = 1, NADJ, 2
          j_neibor = neibor(jb,j)
          ! disconnected graph
          k = D(i,j)
          tmp1 = -up0t_eq(i_neibor,i) - dn0t_eq(i_neibor,i)
          tmp2 =    Gt0up(j,j_neibor) +   Gt0dn(j,j_neibor)
          value_kk(k) = value_kk(k) + tmp1*tmp2
          tmp1 = -up0t_eq(i_neibor,i) - dn0t_eq(i_neibor,i)
          tmp2 =    Gt0up(j_neibor,j) +   Gt0dn(j_neibor,j)
          value_kk(k) = value_kk(k) + tmp1*tmp2
          tmp1 = -up0t_eq(i,i_neibor) - dn0t_eq(i,i_neibor)
          tmp2 =    Gt0up(j,j_neibor) +   Gt0dn(j,j_neibor)
          value_kk(k) = value_kk(k) + tmp1*tmp2
          tmp1 = -up0t_eq(i,i_neibor) - dn0t_eq(i,i_neibor)
          tmp2 =    Gt0up(j_neibor,j) +   Gt0dn(j_neibor,j)
          value_kk(k) = value_kk(k) + tmp1*tmp2
          ! connected graph
          value_kk(k) = value_kk(k) &
                      - upt0_un(j_neibor,i)*up0t_un(i_neibor,j)&
                      - dnt0_un(j_neibor,i)*dn0t_un(i_neibor,j)
          value_kk(k) = value_kk(k) &
                      - upt0_un(j,i)*up0t_un(i_neibor,j_neibor)&
                      - dnt0_un(j,i)*dn0t_un(i_neibor,j_neibor)
          value_kk(k) = value_kk(k) &
                      - upt0_un(j_neibor,i_neibor)*up0t_un(i,j)&
                      - dnt0_un(j_neibor,i_neibor)*dn0t_un(i,j)
          value_kk(k) = value_kk(k) &
                      - upt0_un(j,i_neibor)*up0t_un(i,j_neibor)&
                      - dnt0_un(j,i_neibor)*dn0t_un(i,j_neibor)
       enddo
       enddo
    enddo
    enddo

    ! <K(T_i)D(0)>
    do i  = 1, n    
    do ib = 1, NADJ, 2
       i_neibor = neibor(ib,i)
       do j = 1, n
          k = D(i,j)
          value_kd(k) = value_kd(k) - up0t_eq(i_neibor,i)*Gt0up(j,j)*Gt0dn(j,j)
          value_kd(k) = value_kd(k) - dn0t_eq(i_neibor,i)*Gt0dn(j,j)*Gt0up(j,j)
          value_kd(k) = value_kd(k) - upt0_un(j,i)*up0t_un(i_neibor,j)*Gt0dn(j,j)
          value_kd(k) = value_kd(k) - dnt0_un(j,i)*dn0t_un(i_neibor,j)*Gt0up(j,j)

          value_kd(k) = value_kd(k) - up0t_eq(i,i_neibor)*Gt0up(j,j)*Gt0dn(j,j)
          value_kd(k) = value_kd(k) - dn0t_eq(i,i_neibor)*Gt0dn(j,j)*Gt0up(j,j)
          value_kd(k) = value_kd(k) - upt0_un(j,i_neibor)*up0t_un(i,j)*Gt0dn(j,j)
          value_kd(k) = value_kd(k) - dnt0_un(j,i_neibor)*dn0t_un(i,j)*Gt0up(j,j)
       enddo
    enddo
    enddo

    !<D(T_i)K(0)> and <D(T_i)D(0)>
    do i = 1, n
       do j  = 1, n
       do jb = 1, NADJ, 2
          j_neibor = neibor(jb,j)
          k = D(i,j)
          value_dk(k) = value_dk(k) + upt0_un(j_neibor,i)*up0t_un(i,j)*dn0t_eq(i,i)
          value_dk(k) = value_dk(k) + dnt0_un(j_neibor,i)*dn0t_un(i,j)*up0t_eq(i,i)
          value_dk(k) = value_dk(k) + up0t_eq(i,i)*Gt0up(j,j_neibor)*dn0t_eq(i,i)
          value_dk(k) = value_dk(k) + dn0t_eq(i,i)*Gt0dn(j,j_neibor)*up0t_eq(i,i)

          value_dk(k) = value_dk(k) + upt0_un(j,i)*up0t_un(i,j_neibor)*dn0t_eq(i,i)
          value_dk(k) = value_dk(k) + dnt0_un(j,i)*dn0t_un(i,j_neibor)*up0t_eq(i,i)
          value_dk(k) = value_dk(k) + up0t_eq(i,i)*Gt0up(j_neibor,j)*dn0t_eq(i,i)
          value_dk(k) = value_dk(k) + dn0t_eq(i,i)*Gt0dn(j_neibor,j)*up0t_eq(i,i)
       enddo
       !<D(T_i)D(0)>
       tmp1 = -up0t_eq(i,i)*Gt0up(j,j)-upt0_un(j,i)*up0t_un(i,j)
       tmp2 = -dn0t_eq(i,i)*Gt0dn(j,j)-dnt0_un(j,i)*dn0t_un(i,j)
       value_dd(k) = value_dd(k) + tmp1*tmp2
       enddo
    enddo
 
    !write(*,*) 'In un-equal time subroutine:'
    !write(*,*) 'alphax=',alphax
    !write(*,*) 'KK=',value_kk, ' KD=',value_kd
    !write(*,*) 'DK=',value_dk, ' DD=',value_dd
    !pause
    dynmod_tkk(:) = dynmod_tkk(:) + value_kk(:)
    dynmod_tkd(:) = dynmod_tkd(:) + value_kd(:) + value_dk(:)
    dynmod_tdd(:) = dynmod_tdd(:) + value_dd(:)
    dynmod_t(:)   = dynmod_t(:) + value_kk(:) - alphax*value_kd(:)& 
                    - alphax* value_dk(:) + alphax**2*value_dd(:)

    deallocate(value_kk,value_kd,value_dk,value_dd)

  end subroutine DQMC_TDMX_MeasALL
  !--------------------------------------------------------------------!
  subroutine DQMC_TDMX_EQT(TX,upt0,up0t,dnt0,dn0t,Gt0up,Gt0dn)
    !
    ! Purpose
    ! =======
    ! Measure suscept, <T H_K(t) H_K(0)> for t=0 case
    !
    ! Arguements
    ! ==========
    type(TDMX_RT), intent(inout)    :: TX
    real(wp),   intent(in)       :: upt0(:,:)
    real(wp),   intent(in)       :: up0t(:,:)
    real(wp),   intent(in)       :: dnt0(:,:)
    real(wp),   intent(in)       :: dn0t(:,:)
    real(wp),   intent(inout)    :: Gt0up(:,:)
    real(wp),   intent(inout)    :: Gt0dn(:,:)
    
    ! ... Local scalar ...
    integer  :: n, NADJ, nClass
    integer  :: i,ib, i_neibor
    integer  :: j,jb, j_neibor,k
    ! neighborhood array
    integer,  pointer :: neibor(:,:)
    integer,  pointer :: D(:,:)
    real(wp) :: alphax,tmp1,tmp2,sum_t
    real(wp), allocatable :: value_kk(:), value_kd(:)
    real(wp), allocatable :: value_dk(:), value_dd(:)
    real(wp), pointer :: dynmod_tkk(:),dynmod_tkd(:)
    real(wp), pointer :: dynmod_tdd(:),dynmod_t(:)
    ! ... Executable ...

    ! Initialization
    n        =  TX%n
    nClass   =  TX%nClass
    NADJ     =  TX%NADJ
    D        => TX%D
    neibor   => TX%neibor

    allocate(value_kk(nClass))
    allocate(value_kd(nClass))
    allocate(value_dk(nClass))
    allocate(value_dd(nClass))

    alphax   =  TX%alphax
    dynmod_tkk  => TX%dynmod_tkk(:,1,TX%tmp)
    dynmod_tkd  => TX%dynmod_tkd(:,1,TX%tmp)
    dynmod_tdd  => TX%dynmod_tdd(:,1,TX%tmp)
    dynmod_t    => TX%dynmod_t(:,1,TX%tmp)

    value_kk = ZERO
    value_kd = ZERO
    value_dk = ZERO
    value_dd = ZERO

    !write(*,*) '-------------------------------------------------------'
    !write(*,*) 'In T=0 case'
    !write(*,*) 'upt0(1,2)=',upt0(1,2),' dnt0(1,2)=',dnt0(1,2)
    !write(*,*) 'up0t(1,2)=',up0t(1,2),' dn0t(1,2)=',dn0t(1,2)
    !write(*,*) 'upt0(2,1)=',upt0(2,1),' dnt0(2,1)=',dnt0(2,1)
    !write(*,*) 'up0t(2,1)=',up0t(2,1),' dn0t(2,1)=',dn0t(2,1)

    ! <K(0)K(0)>
    do i  = 1, n
    do ib = 1, NADJ, 2
       i_neibor = neibor(ib,i)
       do j  = 1, n
       do jb = 1, NADJ, 2
          j_neibor = neibor(jb,j)
          k = D(i,j)
          ! disconnected graph
          tmp1 = -up0t(i_neibor,i) - dn0t(i_neibor,i)
          tmp2 = -up0t(j_neibor,j) - dn0t(j_neibor,j)
          value_kk(k) = value_kk(k) + tmp1*tmp2
          tmp1 = -up0t(i_neibor,i) - dn0t(i_neibor,i)
          tmp2 = -up0t(j,j_neibor) - dn0t(j,j_neibor)
          value_kk(k) = value_kk(k) + tmp1*tmp2
          tmp1 = -up0t(i,i_neibor) - dn0t(i,i_neibor)
          tmp2 = -up0t(j_neibor,j) - dn0t(j_neibor,j)
          value_kk(k) = value_kk(k) + tmp1*tmp2
          tmp1 = -up0t(i,i_neibor) - dn0t(i,i_neibor)
          tmp2 = -up0t(j,j_neibor) - dn0t(j,j_neibor)
          value_kk(k) = value_kk(k) + tmp1*tmp2
          ! connected graph
          value_kk(k) = value_kk(k) &
                      - upt0(j_neibor,i)*up0t(i_neibor,j)&
                      - dnt0(j_neibor,i)*dn0t(i_neibor,j)
          value_kk(k) = value_kk(k) &
                      - upt0(j,i)*up0t(i_neibor,j_neibor)&
                      - dnt0(j,i)*dn0t(i_neibor,j_neibor)
          value_kk(k) = value_kk(k) &
                      - upt0(j_neibor,i_neibor)*up0t(i,j)&
                      - dnt0(j_neibor,i_neibor)*dn0t(i,j)
          value_kk(k) = value_kk(k) &
                      - upt0(j,i_neibor)*up0t(i,j_neibor)&
                      - dnt0(j,i_neibor)*dn0t(i,j_neibor)
       enddo
       enddo
    enddo
    enddo

    !<K(0)D(0)>
    do i  = 1, n
    do ib = 1, NADJ, 2
       i_neibor = neibor(ib,i)
       do j = 1, n
          k = D(i,j)
          value_kd(k) = value_kd(k) - up0t(i_neibor,i)*up0t(j,j)*dn0t(j,j)
          value_kd(k) = value_kd(k) - dn0t(i_neibor,i)*dn0t(j,j)*up0t(j,j)
          value_kd(k) = value_kd(k) + upt0(j,i)*up0t(i_neibor,j)*dn0t(j,j)
          value_kd(k) = value_kd(k) + dnt0(j,i)*dn0t(i_neibor,j)*up0t(j,j)

          value_kd(k) = value_kd(k) - up0t(i,i_neibor)*up0t(j,j)*dn0t(j,j)
          value_kd(k) = value_kd(k) - dn0t(i,i_neibor)*dn0t(j,j)*up0t(j,j)
          value_kd(k) = value_kd(k) + upt0(j,i_neibor)*up0t(i,j)*dn0t(j,j)
          value_kd(k) = value_kd(k) + dnt0(j,i_neibor)*dn0t(i,j)*up0t(j,j)
       enddo
    enddo
    enddo



    !<D(0)K(0)> and <D(0)D(0)>
    do i = 1, n
       do j  = 1, n
       do jb = 1, NADJ, 2
          j_neibor = neibor(jb,j)
          k = D(i,j)
          value_dk(k) = value_dk(k) - up0t(i,i)*up0t(j_neibor,j)*dn0t(i,i)
          value_dk(k) = value_dk(k) - dn0t(i,i)*dn0t(j_neibor,j)*up0t(i,i)
          value_dk(k) = value_dk(k) + upt0(j_neibor,i)*up0t(i,j)*dn0t(i,i)
          value_dk(k) = value_dk(k) + dnt0(j_neibor,i)*dn0t(i,j)*up0t(i,i)

          value_dk(k) = value_dk(k) - up0t(i,i)*up0t(j,j_neibor)*dn0t(i,i)
          value_dk(k) = value_dk(k) - dn0t(i,i)*dn0t(j,j_neibor)*up0t(i,i)
          value_dk(k) = value_dk(k) + upt0(j,i)*up0t(i,j_neibor)*dn0t(i,i)
          value_dk(k) = value_dk(k) + dnt0(j,i)*dn0t(i,j_neibor)*up0t(i,i)
       enddo          
       !<D(0)D(0)>
       k = D(i,j)
       tmp1 = up0t(i,i)*up0t(j,j)-upt0(j,i)*up0t(i,j)
       tmp2 = dn0t(i,i)*dn0t(j,j)-dnt0(j,i)*dn0t(i,j)
       value_dd(k) = value_dd(k) + tmp1*tmp2
       enddo
    enddo

    dynmod_tkk(:) = dynmod_tkk(:) + value_kk(:)
    dynmod_tkd(:) = dynmod_tkd(:) + value_kd(:) + value_dk(:)
    dynmod_tdd(:) = dynmod_tdd(:) + value_dd(:)
    dynmod_t(:)   = dynmod_t(:) + value_kk(:) - alphax*value_kd(:)& 
                    - alphax* value_dk(:) + alphax**2*value_dd(:)
    !write(*,*) 'In equal time subroutine:'
    !write(*,*) 'alphax=',alphax
    !write(*,*) 'KK=',value_kk, ' KD=',value_kd
    !write(*,*) 'DK=',value_dk, ' DD=',value_dd

    ! prepare G_kk(0,0) and G_kl for unequal-time susceptibility
    tmp1 = ZERO
    tmp2 = ZERO
    Gt0up = ZERO
    Gt0dn = ZERO
    do i = 1, n
       do ib = 1, NADJ
          i_neibor = neibor(ib,i)
          Gt0up(i,i_neibor) = - up0t(i_neibor,i)
          Gt0dn(i,i_neibor) = - dn0t(i_neibor,i)
          tmp2 = tmp2 - up0t(i_neibor,i) - dn0t(i_neibor,i)
       end do
       Gt0up(i,i) = - up0t(i,i)
       Gt0dn(i,i) = - dn0t(i,i)
       tmp1 = tmp1 + up0t(i,i)*dn0t(i,i)
    end do
    
    TX%dynmod_e(1,TX%tmp) = TX%dynmod_e(1,TX%tmp) + tmp2
    TX%dynmod_e(2,TX%tmp) = TX%dynmod_e(2,TX%tmp) + tmp1


    deallocate(value_kk,value_kd,value_dk,value_dd)

  end subroutine DQMC_TDMX_EQT
  !--------------------------------------------------------------------!
  subroutine DQMC_TDMX_FTM(dtau,TX)
    !
    ! Purpose
    ! =======
    !    This subroutine obtain X(T) and Fourier transform to X(iw_n)
    !
    ! Arguments
    ! =========
    !
    type(TDMX_RT), intent(inout)    :: TX
    real(wp), intent(in)         :: dtau

    integer :: i,j,n,L,nClass
    integer, pointer  :: F(:)
    real(wp) :: tau_x,hx,factor

    n = TX%n
    L = TX%L
    nClass= TX%nClass
    F     => TX%F

    ! Calculate X(\tau)
    do j = 1, L
    do i = 1, nClass
       factor = 1.d0/DBLE(F(i))
       TX%gnl(i,j,TX%tmp) = TX%gnl(i,j,TX%tmp)*factor
       TX%dynmod_tkk(i,j,TX%tmp)=TX%dynmod_tkk(i,j,TX%tmp)*factor
       TX%dynmod_tkd(i,j,TX%tmp)=TX%dynmod_tkd(i,j,TX%tmp)*factor
       TX%dynmod_tdd(i,j,TX%tmp)=TX%dynmod_tdd(i,j,TX%tmp)*factor
       TX%dynmod_t(i,j,TX%tmp)=TX%dynmod_t(i,j,TX%tmp)*factor
    enddo
    enddo
    TX%dynmod_e(:,TX%tmp)=TX%dynmod_e(:,TX%tmp)/DBLE(n)
    ! Calculate X(T)
    TX%dynmod_w0(:,TX%tmp)=ZERO
    hx=2.d0/3.d0
    do i = 1, L
       TX%dynmod_w0(1,TX%tmp)=TX%dynmod_w0(1,TX%tmp) + &
    hx*SUM(TX%dynmod_tkk(:,i,TX%tmp))
       TX%dynmod_w0(2,TX%tmp)=TX%dynmod_w0(2,TX%tmp) + &
    hx*SUM(TX%dynmod_tkd(:,i,TX%tmp))
       TX%dynmod_w0(3,TX%tmp)=TX%dynmod_w0(3,TX%tmp) + &
    hx*SUM(TX%dynmod_tdd(:,i,TX%tmp))
       TX%dynmod_w0(4,TX%tmp)=TX%dynmod_w0(4,TX%tmp) + &
    hx*SUM(TX%dynmod_t(:,i,TX%tmp))
            hx=2.d0-hx
    end do
    TX%dynmod_w0(:,TX%tmp)=TX%dynmod_w0(:,TX%tmp)*dtau

  end subroutine DQMC_TDMX_FTM
  !--------------------------------------------------------------------!

  subroutine DQMC_TDMX_Acc(TX,sgn)
    !
    ! Purpose
    ! =======
    !    This subroutine accumulates TDMX.
    !
    ! Arguments
    ! =========
    !
    type(TDMX_RT), intent(inout) :: TX
    real(wp),   intent(in)    :: sgn

    ! ... Local scalar ...

    integer  :: i, j, L
    real(wp), pointer :: dynmod_tkk(:,:),dynmod_tkd(:,:)
    real(wp), pointer :: dynmod_t(:,:),dynmod_tdd(:,:)
    real(wp), pointer :: dynmod_w0(:)
    real(wp), pointer :: dynmod_e(:)
    real(wp), pointer :: dynmod_gnl(:,:)

    ! ... Executable ...

    ! Initialization
    L     = TX%L
    
    dynmod_gnl => TX%gnl(:,:,TX%tmp)

    dynmod_tkk => TX%dynmod_tkk(:,:,TX%tmp)
    dynmod_tkd => TX%dynmod_tkd(:,:,TX%tmp)
    dynmod_tdd => TX%dynmod_tdd(:,:,TX%tmp)
    dynmod_t => TX%dynmod_t(:,:,TX%tmp)

    dynmod_w0  => TX%dynmod_w0(:,TX%tmp)
    dynmod_e   => TX%dynmod_e(:,TX%tmp)

    ! Accumulate
    TX%gnl(:,:,TX%avg) = TX%gnl(:,:,TX%avg) + dynmod_gnl(:,:)*sgn
    TX%dynmod_tkk(:,:,TX%avg) = TX%dynmod_tkk(:,:,TX%avg) + dynmod_tkk(:,:)*sgn
    TX%dynmod_tkd(:,:,TX%avg) = TX%dynmod_tkd(:,:,TX%avg) + dynmod_tkd(:,:)*sgn
    TX%dynmod_tdd(:,:,TX%avg) = TX%dynmod_tdd(:,:,TX%avg) + dynmod_tdd(:,:)*sgn
    TX%dynmod_t(:,:,TX%avg) = TX%dynmod_t(:,:,TX%avg) + dynmod_t(:,:)*sgn

    TX%dynmod_w0(:,TX%avg) = TX%dynmod_w0(:,TX%avg) + dynmod_w0(:)*sgn
    TX%dynmod_e(:,TX%avg)  = TX%dynmod_e(:,TX%avg)  + dynmod_e(:)*sgn

    !write(*,*) 'BOW=', dynmod_e

    ! Clean up for next measure
    TX%gnl(:,:,TX%tmp) = ZERO
    TX%dynmod_tkk(:,:,TX%tmp) = ZERO
    TX%dynmod_tkd(:,:,TX%tmp) = ZERO
    TX%dynmod_tdd(:,:,TX%tmp) = ZERO
    TX%dynmod_t(:,:,TX%tmp) = ZERO

    TX%dynmod_w0(:,TX%tmp) = ZERO
    TX%dynmod_e(:,TX%tmp)  = ZERO

  end subroutine DQMC_TDMX_Acc

  !--------------------------------------------------------------------!

  subroutine DQMC_TDMX_Avg(TX)
    !
    ! Purpose
    ! =======
    !    This subroutine averges the pair measurements.
    !
    ! Arguments
    ! =========
    !
    type(TDMX_RT), intent(inout) :: TX ! TX

    ! ... BLAS function ...
    real(wp), external :: ddot

    ! ... local scalar ...
    real(wp) :: factor
    real(wp) :: alphax
    integer  :: nl
    ! temperary variable
    integer :: i,L

    ! ... Executable ...
    factor =  ONE/DBLE(TX%cnt)
    nl =  TX%L*TX%nClass
    
    call dscal(nl, factor, TX%gnl(:,:,TX%avg),1)
    call dscal(nl, factor, TX%dynmod_tkk(:,:,TX%avg), 1)
    call dscal(nl, factor, TX%dynmod_tkd(:,:,TX%avg), 1)
    call dscal(nl, factor, TX%dynmod_tdd(:,:,TX%avg), 1)
    call dscal(nl, factor, TX%dynmod_t(:,:,TX%avg),   1)

    TX%dynmod_w0(:,TX%avg) = TX%dynmod_w0(:,TX%avg)*factor
    TX%dynmod_e(:,TX%avg)  = TX%dynmod_e(:,TX%avg)*factor

    !TX%dynmod_tkk(:,:,TX%avg) = TX%dynmod_tkk(:,:,TX%avg)& 
    ! -(TX%dynmod_e(1,TX%avg))**2
    !TX%dynmod_tkd(:,:,TX%avg) = TX%dynmod_tkd(:,:,TX%avg)&
    ! - 2.d0*TX%dynmod_e(1,TX%avg)*TX%dynmod_e(2,TX%avg)
    !TX%dynmod_tdd(:,:,TX%avg) = TX%dynmod_tdd(:,:,TX%avg)&
    ! -(TX%dynmod_e(2,TX%avg))**2
    !TX%dynmod_t(:,:,TX%avg) = TX%dynmod_t(:,:,TX%avg)&
    ! -(TX%dynmod_e(1,TX%avg))**2&
    ! +2.d0*TX%dynmod_e(1,TX%avg)*TX%dynmod_e(2,TX%avg)&
    ! -(TX%dynmod_e(2,TX%avg))**2
    
    !write(*,*) 'The average of per site is:'
    !write(*,*) '<KK>=', TX%dynmod_tkk(:,TX%avg)
    !write(*,*) '<KD>=', TX%dynmod_tkd(:,TX%avg)
    !write(*,*) '<DD>=', TX%dynmod_tdd(:,TX%avg)
    !write(*,*) '<ALL>=', TX%dynmod_t(:,TX%avg)
    !write(*,*) '<K>=', TX%dynmod_e(1,TX%avg)
    !write(*,*) '<D>=', TX%dynmod_e(2,TX%avg)
    !pause

    !For sign
    TX%sgn = TX%sgn*factor

    ! reset counts
    TX%cnt = 0

  end subroutine DQMC_TDMX_Avg

  !--------------------------------------------------------------------!

  subroutine DQMC_TDMX_WriteBin(TX,idx,rank)
  !
  ! Purpose
  ! =======
  !    This subroutine averges the pair measurements.
  !
  ! Arguments
  ! =========
  ! rank: processor index for parallel job, single processor rank=0
  !
  type(TDMX_RT), intent(inout) :: TX ! TX
  integer, intent(in)    :: rank
  integer, intent(in)    :: idx
  ! temperary variable
  integer :: ic,j,L,Nc
  integer :: unit_t
  real(wp)::alphax
 
  L = TX%L
  Nc = TX%nClass
  alphax = TX%alphax

  !do j=1, L
  !   TX%dynmod_tkk(1,j,TX%avg)=sum(TX%dynmod_tkk(:,j,TX%avg))
  !enddo
  unit_t = 10000+rank
  write(unit_t,"('#No. bin=',i5,' from processor',i5)") idx,rank
  write(unit_t,"(6(1X,f12.8))") ((TX%dynmod_tkk(ic,j,TX%avg),j=1,L),ic=1,Nc)&
  ,TX%sgn,TX%dynmod_e(1,TX%avg)

  unit_t = 20000+rank
  write(unit_t,"('#No. bin=',i5,' from processor',i5)") idx,rank
  write(unit_t,"(6(1X,f12.8))") ((TX%dynmod_tkd(ic,j,TX%avg),j=1,L),ic=1,Nc)&
  ,TX%sgn,DSQRT(TX%dynmod_e(1,TX%avg)*TX%dynmod_e(2,TX%avg)*2.d0)

  !do j=1, L
  !   TX%dynmod_tdd(1,j,TX%avg)=sum(TX%dynmod_tdd(:,j,TX%avg))
  !enddo
  unit_t = 30000+rank
  write(unit_t,"('#No. bin=',i5,' from processor',i5)") idx,rank
  write(unit_t,"(6(1X,f12.8))") ((TX%dynmod_tdd(ic,j,TX%avg),j=1,L),ic=1,Nc)&
  ,TX%sgn,TX%dynmod_e(2,TX%avg)

  !do j=1, L
  !   TX%dynmod_t(1,j,TX%avg)=sum(TX%dynmod_t(:,j,TX%avg))
  !enddo
  unit_t = 40000+rank
  write(unit_t,"('#No. bin=',i5,' from processor',i5)") idx,rank
  write(unit_t,"(6(1X,f12.8))") ((TX%dynmod_t(ic,j,TX%avg),j=1,L),ic=1,Nc)&
  ,TX%sgn,-TX%dynmod_e(1,TX%avg)+TX%dynmod_e(2,TX%avg)*alphax

  unit_t = 50000+rank
  write(unit_t,"('#No. bin=',i5,' from processor',i5)") idx,rank
  write(unit_t,"(6(1X,f12.8))") ((TX%gnl(ic,j,TX%avg),j=1,L),ic=1,Nc)&
  ,TX%sgn, 0.d0

  TX%gnl = ZERO
  TX%dynmod_tkk = ZERO
  TX%dynmod_tkd = ZERO
  TX%dynmod_tdd = ZERO
  TX%dynmod_t = ZERO

  TX%dynmod_w0 = ZERO
  TX%dynmod_e = ZERO
  TX%sgn = ZERO

  end subroutine DQMC_TDMX_WRITEBIN
  !--------------------------------------------------------------------!

  subroutine DQMC_Make_FTMX(L, itvl, dtau, TX, base)
    !
    ! Purpose
    ! =======
    !    This subroutine constructs Fourier transformation matrix (FTM).
    !    It implements numerical itegration on FT.
    !
    ! Arguments
    ! =========
    !
    type(TDMX_RT), intent(inout)   :: TX
    real(wp), intent(in)      :: base, dtau
    integer, intent(in)       :: L
    integer, intent(in)       :: itvl

    ! ... parameters ...
    real(wp), parameter  :: ONETHIRD  = ONE/3.0_wp
    real(wp), parameter  :: TWOTHIRD  = TWO/3.0_wp
    real(wp), parameter  :: FOURTHIRD = 4.0_wp/3.0_wp
    real(wp), parameter  :: TWOPI     = 6.283185307179586
    complex(wp), parameter  :: EI     = cmplx(ZERO, ONE)

    ! ... Local variables ...
    real(wp),pointer :: FTM(:,:)
    integer  :: i, j, img_idx, longL
    real(wp) :: r
    complex(wp) :: coef, omega


    ! ... Executable ...
    FTM =>TX%FTM
    ! r is just dt
    r       = dtau / itvl
    img_idx = L
    longL   = itvl*(L-1)+1

    do j = 1, L
       ! The first term 
       FTM(1, j) = ONETHIRD * r
       FTM(1, j+img_idx) = ZERO

       ! The middle terms
       omega = EI*(j-base)*TWOPI/L/itvl
       do i = 2, longL-1, 2
          coef = FOURTHIRD * exp(omega*(i-1)) * r
          FTM(i,j)         = real (coef)          ! real part
          FTM(i,j+img_idx) = aimag(coef)          ! imagianry part

          coef =  TWOTHIRD * exp(omega*i) * r
          FTM(i+1,j)         = real (coef)        ! real part
          FTM(i+1,j+img_idx) = aimag(coef)        ! imagianry part
       end do

       ! The last second term
       if (mod(longL-1,2) .ne. 0) then
          coef = FOURTHIRD * exp(omega*(longL-2)) * r
          FTM(longL-1, j)        = real (coef)    ! real part
          FTM(longL-1,j+img_idx) = aimag(coef)    ! imagianry part
       end if

       ! The last term
       coef = ONETHIRD * exp(omega*(longL-1)) * r
       FTM(longL,j)         = real (coef)         ! real part
       FTM(longL,j+img_idx) = aimag(coef)         ! imagianry part
    end do

  end subroutine DQMC_Make_FTMX

  !--------------------------------------------------------------------!
  subroutine DQMC_SplineX(L, dtau, itvl, input, output, W)
    !  
    ! Purpose
    ! =======
    !    This subroutine computes spline interpolation of a set of data.
    !    Fx is a nxL matrix where L is the time slice. For the details
    !    of spline interpolation, see [1]
    !
    !    1. Solve a tridiagonal linear system (L-2)x(L-2)
    !       The following fomula only shows the first row of Fx. 
    !       Let F(i) = Fx(1,i)
    !       
    !                 [4 1       ]          [   F(1)-2F(2)+F(3)  ]
    !                 [1 4 1     ]      6   [   F(2)-2F(3)+F(4)  ]
    !               X [  ......  ] = ------ [          :         ]
    !                 [     1 4 1]   dtau^2 [          :         ]
    !                 [       1 4]          [ F(L-2)-2F(L-1)+F(L)]
    !
    !    2. Formulate a_i*dx^3 + b_i*dx^2+c_i*dx+d_i
    !       where a_i = (X_{i+1}-M_i)/(6*dtau)
    !             b_i = X_i/2
    !             c_i = (F(i+1)-F(i))/dtau - (X_{i+1}+2X_{i})*dtau/6
    !             d_i = F(i)
    !
    !   [1] Van Loan, Charles F. Introduction to Scientific Computing. 
    !       1997 New Jersey: Prentice Hall
    !
    ! Arguments
    ! =========
    !
    integer,  intent(in)    :: L
    real(wp), intent(in)    :: dtau
    real(wp), intent(in)    :: input(:)         ! original data
    real(wp), target, intent(inout) &
                      :: output(itvl*(L-1)+1)   ! interpolated points
    real(wp), intent(inout) :: W(:)               ! working space
    integer,  intent(in)    :: itvl               ! itvl must >= 2

    ! ... Local ...
    real(wp) :: r, h               ! temp variable
    integer  :: i, j, k
    real(wp) :: a, b, c, d         ! coefficient of spline
    real(wp), pointer :: Pt(:)

    ! ... Executable ...

    ! We use the Thomas algorithm to solve the tridiagonal linear
    ! system.

    !! Reuse the last block of the output array as a working space
    Pt => output((itvl-1)*(L-1)+1:itvl*(L-1)+1)

    r = 6.0_wp/dtau/dtau
    do i = 1,L-2
       Pt(i+1) = r*(input(i)-TWO*input(i+1)+input(i+2))
    end do

    ! forward subsitution
    W(2) = 4.0_wp
    do i = 3,L-1
       r     = ONE/W(i-1)
       W(i)  = 4.0_wp - r
       Pt(i) = Pt(i) - r*Pt(i-1)
    end do

    ! back subsitution
    r = ONE/W(L-1)
    Pt(1)   = ZERO
    Pt(L)   = ZERO
    Pt(L-1)  = r*Pt(L-1)

    do i = L-2,2,-1
       r     = ONE/W(i)
       Pt(i) = r*(Pt(i)-Pt(i+1))
    end do

    ! Find interpolated points

    ! There are L-1 curves for L points
    r = dtau/itvl        ! step size
    do i = 1, L-1
       a = (Pt(i+1)-Pt(i))/6.0_wp/dtau
       b = Pt(i)/TWO
       c = (input(i+1)-input(i))/dtau-&
               (Pt(i+1)+2*Pt(i))*dtau/6.0_wp
       d = input(i)

       ! The Ft(itvl*(i-1)+1) = input(i) 
       output(itvl*(i-1)+1) = input(i)
       h = r
       do k = itvl*(i-1)+2, itvl*i
          output(k) = ((a*h+b)*h+c)*h+d
          h = h + r
       end do
    end do

    ! The Ft(itvl*L) = Pt(L) 
    output(itvl*(L-1)+1) = input(L)

  end subroutine DQMC_SplineX

  !--------------------------------------------------------------------!


  subroutine DQMC_FT_TimeX(L, itvl, dtau, data_t, data_w, FTM)
    !  
    ! Purpose
    ! =======
    !    This subroutine computes FT on given data.
    !  
    !   
    ! Arguments
    ! =========
    !
    integer,  intent(in)    :: L             ! number of time slices
    real(wp), intent(in)    :: data_t(:)     ! L+2 input, time
    real(wp), intent(inout) :: data_w(:)     ! L+2 ourput,frequency
    real(wp), intent(in)    :: FTM(:, :)     ! longL*(L+2), Fourier matrix
    real(wp), intent(in)    :: dtau          ! step size
    integer,  intent(in)    :: itvl          ! spline interval

    ! ... local varaibles ...
    real(wp),pointer :: data_tlong(:) ! longL, output of spline
    real(wp),pointer :: W(:)          ! working space: dim=L+1
    integer  :: i,longL

    ! ... Executable ...
    longL    = (L-1)*itvl+1
    !allocate work array
    allocate(data_tlong(longL))
    allocate(W(L+1))

    ! Spine dynmod_t
    call DQMC_SplineX(L, dtau, itvl, data_t, data_tlong, W)
    !do i = 1, longL
    !   write(*,*) i, data_tlong(i)
    !enddo
    ! FT
    call dgemm('N', 'N', 1, 2*L, longL, ONE, data_tlong, 1, &
                FTM, longL, ZERO, data_w, 1)

    !do i = 1, L
    !   write(*,*) i, data_w(i)
    !enddo
    !pause
    deallocate(data_tlong,W)

  end subroutine DQMC_FT_TimeX

end module DQMC_TDMX_RT
