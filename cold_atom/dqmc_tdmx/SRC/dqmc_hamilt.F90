module DQMC_HAMILT

  use DQMC_GEOM_PARAM
  use DQMC_LATT
  use DQMC_RECLATT

  implicit none

  type :: Hamiltonian_t
     integer                :: nJsites,nUsites,ntsites            !number of sites having a non-zero J,U,t
     integer                :: maxneig                            !maximum number of neighbors
     integer, pointer       :: tnneig(:),Unneig(:),Jnneig(:)      !Number of neighbors of each site (0:nsites-1)
     integer, pointer       :: tsite(:),Usite(:),Jsite(:)         !Sites havind a non-zero J,U,t (0:nsites-1)
     integer, pointer       :: tneig(:,:),Uneig(:,:),Jneig(:,:)   !Sites neighboring each site   (0:nsites-1,nsites)

     real*8                 :: mu                                 !chemical potential
     real*8, pointer        :: Uv(:,:),Jv(:,:)                    !value of U and J for each pair (0:nsites-1, 0:nsites-1)
     real*8, pointer        :: Uvalue(:),muvalue(:)               !values of U and mu inequivalent by symmetry (nlocclass)

     complex*16, pointer    :: hop(:,:)                           !value of t for each pair (0:nsites-1, 0:nsites-1)
     complex*16, pointer    :: tvalue(:)                          !values of t inequivalent by symmetry (nhopclass) 
     complex*16, pointer    :: phase(:)                           !wave function phase(not used in QMC)


     integer                :: nhopclass                          !number of different hoppings
     integer                :: nlocclass                          !number of sites with different U/mu
     integer, pointer       :: mylocclass(:)                      !class for each site having U/mu (0:nsites-1)
     integer, pointer       :: myhopclass(:,:)                    !hopping class for each pair of neighboring site (0:nsites-1,maxval(tnneig))

     logical                :: constructed
     logical                :: neig_found
     logical                :: analyzed
  end type Hamiltonian_t

contains

  !--------------------------------------------------------------------------------------
  ! Fill the hamiltonian
  !--------------------------------------------------------------------------------------
  subroutine construct_hamilt(hamilt,lattice,recip_lattice)
    type(lattice_t),intent(in)       :: lattice
    type(recip_lattice_t),intent(in) :: recip_lattice
    type(Hamiltonian_t),intent(out)    :: hamilt
    integer              :: maxhop, iat,jat,natom,nsites,ios,icount,ihop,j
    real*8               :: hop3d(rdim), tijtmp, Jtmp, Utmp,ktwist(rdim), kpoint(rdim)
    real*8, pointer      :: Uv(:,:),Jv(:,:)
    complex*16, pointer  :: hop(:,:),phase(:)
    integer, allocatable :: nhop(:),hopto(:,:),hoptarget(:,:)
    real*8, allocatable  :: tij(:,:),Uat(:,:),Jint(:,:),twisthop(:,:),ahop(:,:,:)
    character*50         :: string
    logical              :: ldum


    if(.not.lattice%constructed)stop'Need to construct lattice before building Hamiltonian'
    if(.not.recip_lattice%initialized)stop'Need to initialize recip_lattice before building Hamiltonian'

    natom=lattice%natom
    nsites=lattice%nsites
    ktwist(1:rdim) = recip_lattice%ktwist(1:rdim)
    kpoint(1:rdim) = recip_lattice%kpoint(1:rdim)
    allocate(hop(0:nsites-1,0:nsites-1),Uv(0:nsites-1,0:nsites-1),Jv(0:nsites-1,0:nsites-1),phase(0:nsites-1))

    !Read chemical potential
    ldum=move_to_record(INPUT_FIELDS(MU_F),inpunit)
    read(inpunit,*,iostat=j)hamilt%mu
    if(j/=0)stop 'Problem reading #MU'

    !Read the hamiltonian
    ldum=move_to_record(INPUT_FIELDS(HAMILT_F),inpunit)
    allocate(nhop(0:natom-1))
    call count_nhop(nhop,natom); maxhop=maxval(nhop)
    allocate(hopto(0:nsites-1,maxhop),hoptarget(0:natom-1,maxhop),           &
         & tij(0:natom-1,maxhop),Uat(0:natom-1,maxhop),twisthop(0:nsites-1,maxhop),    &
         & Jint(0:natom-1,maxhop),ahop(rdim,0:natom-1,maxhop))
    nhop(:)=0
    ldum=move_to_record(INPUT_FIELDS(HAMILT_F),inpunit)
    do
       read(inpunit,'(A)')string
       read(string,*,iostat=ios)iat,jat,hop3d(1),hop3d(2),hop3d(3),tijtmp,Jtmp,Utmp
       if(ios.ne.0)exit
       nhop(iat)=nhop(iat)+1
       hoptarget(iat,nhop(iat))=jat
       ahop(1,iat,nhop(iat))=hop3d(1)
       ahop(2,iat,nhop(iat))=hop3d(2)
       ahop(3,iat,nhop(iat))=hop3d(3)
       tij(iat,nhop(iat))=tijtmp
       Jint(iat,nhop(iat))=Jtmp
       Uat(iat,nhop(iat))=Utmp
       if(jat.eq.iat.and.sum(hop3d**2)<1.d-10)cycle
       nhop(jat)=nhop(jat)+1
       hoptarget(jat,nhop(jat))=iat
       ahop(1,jat,nhop(jat))=-hop3d(1)
       ahop(2,jat,nhop(jat))=-hop3d(2)
       ahop(3,jat,nhop(jat))=-hop3d(3)
       tij(jat,nhop(jat))=tijtmp
       Jint(jat,nhop(jat))=Jtmp
       Uat(jat,nhop(jat))=Utmp
    enddo

    do iat=0,nsites-1
       jat=mod(iat,natom)
       phase(iat)=exp(im*sum((kpoint(:)-ktwist(:))*(lattice%cartpos(:,iat)-lattice%cartpos(:,jat))))* &
            exp(im*(-sum(ktwist(:)*lattice%cartpos(:,iat))))
    enddo

    !construct hopping table
    do iat=0,nsites-1
       icount=mod(iat,natom)
       do ihop=1,nhop(icount)
          twisthop(iat,ihop)=-sum(ktwist(:)*ahop(:,icount,ihop))
          jat=hoptarget(icount,ihop)
          hopto(iat,ihop)=hoptowho(iat,ahop(1,icount,ihop),jat,lattice)
       enddo
    enddo

    !construct hamiltonian matrix for t,U and J
    hop(:,:)=0.d0; Jv(:,:)=0.d0; Uv(:,:)=0.d0
    do iat=0,nsites-1
       !iat -> neigbors
       icount=mod(iat,natom)
       do ihop=1,nhop(icount)
          j=hopto(iat,ihop)
          !hop is the hopping matrix element <j|T|iat> for c_j+ c_iat
          hop(iat,j)=hop(iat,j)+tij(icount,ihop)*exp(im*twisthop(iat,ihop))
          Jv(iat,j)=Jv(iat,j)+Jint(icount,ihop)
          Uv(iat,j)=Uv(iat,j)+Uat(icount,ihop)
       enddo
    enddo
    hamilt%Uv => Uv
    hamilt%Jv => Jv
    hamilt%hop => hop
    hamilt%phase => phase
    deallocate(nhop,hopto,hoptarget,tij,Uat,twisthop,Jint,ahop)
    hamilt%constructed=.true.

    call find_neighbors(hamilt)
    hamilt%neig_found=.true.

  end subroutine construct_hamilt




  !--------------------------------------------------------------------
  !Count the number of hops that are possible from each orbitals 
  !in the primitive cell
  !--------------------------------------------------------------------
  subroutine count_nhop(nhop,natom)
    integer, intent(in)  :: natom
    integer, intent(out) :: nhop(0:natom-1)
    integer              :: iat,jat,ios
    character*50         :: string
    real*8               :: hop(rdim),tijtmp,Jtmp,Utmp
    do iat=0,natom-1
       nhop(iat)=0
    enddo
    do
       read(inpunit,'(A)')string
       read(string,*,iostat=ios)iat,jat,hop(1),hop(2),hop(3),tijtmp,Jtmp,Utmp
       if(ios.ne.0)exit
       if(iat>natom-1.or.jat>natom-1)then
          write(*,*)'One of the atom in hopping is unspecified'; stop
       endif
       nhop(iat)=nhop(iat)+1
       if(jat.eq.iat.and.sum(hop**2)<1.d-10)cycle
       nhop(jat)=nhop(jat)+1
    enddo
    rewind(inpunit)
  end subroutine count_nhop




  !---------------------------------------------------------------------------------
  ! Given the matrices Uv,Jv and hop it finds which atoms are neighbors and classify 
  !neighbors as with respect to interaction or hopping. Returns tnneig(is): number
  !of neighbors of is; tneig(is,j): j-th neighbor of "is". Analogous definition 
  !holds for Unneig,Uneig and Jnneig,Jneig.
  !------------------------------------------------------------------------------
  subroutine find_neighbors(hamilt)
    integer is,js,ntsites,nUsites,nJsites,nsites
    integer, pointer :: tnneig(:),Unneig(:),Jnneig(:),tneig(:,:),Uneig(:,:),Jneig(:,:),tsite(:),Usite(:),Jsite(:)
    type(Hamiltonian_t),intent(inout) :: hamilt
    !
    if(.not.hamilt%constructed)stop'Hamiltonian needs to be constructed before neig can be found'
    !
    nsites=size(hamilt%hop,1)
    if(associated(tnneig))deallocate(tnneig,Unneig,Jnneig,tneig,Uneig,Jneig,tsite,Usite,Jsite)
    allocate(tnneig(0:nsites-1),Unneig(0:nsites-1),Jnneig(0:nsites-1),  &
         & tneig(0:nsites-1,nsites),Uneig(0:nsites-1,nsites),                 &
         & Jneig(0:nsites-1,nsites),tsite(nsites),Usite(nsites),Jsite(nsites))
    do is=0,nsites-1
       tnneig(is)=0; Unneig(is)=0; Jnneig(is)=0
       do js=0,nsites-1
          if(abs(hamilt%hop(is,js)).gt.1d-9.and.is/=js)then
             tnneig(is)=tnneig(is)+1
             tneig(is,tnneig(is))=js
          endif
          if(abs(hamilt%Uv(is,js)).gt.1d-9)then
             Unneig(is)=Unneig(is)+1
             Uneig(is,Unneig(is))=js
          endif
          if(abs(hamilt%Jv(is,js)).gt.1d-9)then
             Jnneig(is)=Jnneig(is)+1
             Jneig(is,Jnneig(is))=js
          endif
       enddo
    enddo
    ntsites=0; nUsites=0; nJsites=0
    do is=0,nsites-1
       if(tnneig(is).gt.0)then
          ntsites=ntsites+1
          tsite(ntsites)=is
       endif
       if(Unneig(is).gt.0)then
          nUsites=nUsites+1
          Usite(nUsites)=is
       endif
       if(Jnneig(is).gt.0)then
          nJsites=nJsites+1
          Jsite(nJsites)=is
       endif
    enddo
    write(*,*)'Hamiltonian info'
    write(*,*)
    if(maxval(tnneig).eq.0)then
       write(*,*)'Particles are not hopping' 
    else
       write(*,*)'Neighbours tables for hopping.'
       do is=0,nsites-1
          write(*,100)is,(tneig(is,js),hamilt%hop(is,tneig(is,js)),js=1,tnneig(is))
       enddo
    endif
    write(*,*)
    if(maxval(unneig).eq.0)then
       write(*,*)'No density-density interaction.' 
    else
       write(*,*)'Neighbours tables for U'
       do is=0,nsites-1
          write(*,101)is,(Uneig(is,js),hamilt%Uv(is,Uneig(is,js)),js=1,Unneig(is))
       enddo
    endif
    write(*,*)
    if(maxval(jnneig).eq.0)then
       write(*,*)'No spin-spin interaction.' 
    else
       write(*,*)'Neighbours tables for J'
       do is=0,nsites-1
          write(*,101)is,(Jneig(is,js),hamilt%Jv(is,Jneig(is,js)),js=1,Jnneig(is))
       enddo
    endif
    !Fill hamiltonian variables
    hamilt%maxneig=max(maxval(tnneig),maxval(Unneig),maxval(Jnneig))
    !number of neighbors of each site
    hamilt%tnneig=>tnneig
    hamilt%Unneig=>Unneig
    hamilt%Jnneig=>Jnneig
    !neighbors of each site
    hamilt%tneig=>tneig
    hamilt%Uneig=>Uneig
    hamilt%Jneig=>Jneig
    !number of sites on which either t or U or J is different from 0
    hamilt%ntsites=ntsites
    hamilt%nUsites=nUsites
    hamilt%nJsites=nJsites
    !sites on which either t or U or J is different from 0
    hamilt%tsite=>tsite
    hamilt%Usite=>Usite
    hamilt%Jsite=>Jsite
    write(*,*)'===================================================================='
    return
100 format(1x,i4,20(1x,i4,1x,f7.3,1x,f7.3))
101 format(1x,i4,20(1x,i4,1x,f7.3))
  end subroutine find_neighbors

  !----------------------------------------------------------------------------------
  ! Construct a list of only hopping classes. These are a subset of all the
  ! distance classes for which the hopping is non zero.
  !  Find which hoppings are equivalent by symmetry. nhopsite(iat,iclass) returns how many
  ! hoppings in class iclass iat has. ordered_neig(iat,ineig) returns the ineig-th
  ! neighbor of iat. Neighbors are here ordered according to classes. So if
  ! nhopsite(2,1)=3 the first three neigbors of 2 listed in ordered_neig are those
  ! belonging to class 1. tneig has a similar content but neighbors are not ordered.
  !----------------------------------------------------------------------------------
  subroutine count_hop_class(lattice,hamilt)
    type(Hamiltonian_t) :: hamilt
    type(lattice_t) :: lattice
    integer              :: isite,ineig,jsite,newclass,ihc,maxtneig,nsites,hopclass(lattice%nclass),nhopclass
    integer, allocatable :: pairhopclass(:,:)
    complex*16           :: tvaluetmp(lattice%nclass)
    nsites=size(hamilt%tnneig)
    nhopclass=0
    hopclass(:)=0
    allocate(pairhopclass(0:nsites-1,nsites))
    pairhopclass(:,:)=-1
    !Select classes for which hopping is different from 0
    do isite=0,nsites-1
       do ineig=1,hamilt%tnneig(isite)
          jsite=hamilt%tneig(isite,ineig) 
          if(isite==jsite)cycle
          newclass=lattice%myclass(isite,jsite)
          !First check that the class has not been already found
          do ihc=1,nhopclass
             if(newclass==hopclass(ihc))exit
          enddo
          !If hopclass is new, increase number of classes
          if(ihc==nhopclass+1)then 
             hopclass(ihc)=newclass
             nhopclass=ihc
             tvaluetmp(ihc)=hamilt%hop(isite,jsite)
          endif
          !Assign pair to the class
          pairhopclass(isite,ineig)=ihc
       enddo
    enddo
    hamilt%nhopclass=nhopclass
    maxtneig=maxval(hamilt%tnneig(:))
    allocate(hamilt%myhopclass(0:nsites-1,maxtneig),hamilt%tvalue(nhopclass))
    hamilt%tvalue(1:nhopclass)=tvaluetmp(1:nhopclass)
    hamilt%myhopclass(0:nsites-1,1:maxtneig)=pairhopclass(0:nsites-1,1:maxtneig)
    !allocate(nhopsite(0:nsites-1,nhopclass),ordered_neig(0:nsites-1,maxtneig),tvalue(nhopclass))
    !tvalue(1:nhopclass)=tvaluetmp(1:nhopclass)
    !nhopsite(:,:)=0
    !do isite=0,nsites-1
    ! do ineig=1,tnneig(isite)
    !  ihc=pairhopclass(isite,ineig)
    !  nhopsite(isite,ihc)=nhopsite(isite,ihc)+1
    ! enddo
    ! ifound=0
    ! do ihc=1,nhopclass
    !  do ineig=1,tnneig(isite)
    !   if(pairhopclass(isite,ineig)==ihc)then
    !    ifound=ifound+1
    !    ordered_neig(isite,ifound)=tneig(isite,ineig)
    !   endif
    !  enddo
    ! enddo
    !enddo
    !write(*,*)'nhopclass',nhopclass
    !do isite=0,nsites-1
    ! write(*,'(100i4)')isite,(nhopsite(isite,ihc),ihc=1,nhopclass)
    !enddo
    !write(*,*)'Ordered neighbors'
    !do isite=0,nsites-1
    ! write(*,'(100i4)')isite,(ordered_neig(isite,ineig),ineig=1,tnneig(isite))
    !enddo
    deallocate(pairhopclass)
  end subroutine count_hop_class




  !--------------------------------------------------------------------------
  ! Construct a list of only local classes i.e. the subset of the
  ! distance classes defined on the same site. Uvalue contains the U
  ! for that class.
  !--------------------------------------------------------------------------
  subroutine count_local_classes(lattice,hamilt)
    type(Hamiltonian_t) :: hamilt
    type(lattice_t) :: lattice
    integer isite,iclass,localtmp(lattice%nclass),jclass,nsites,nlocclass
    real*8 Utmp(lattice%nclass),mutmp(lattice%nclass)
    nlocclass=0
    nsites=size(hamilt%tnneig)
    allocate(hamilt%mylocclass(0:nsites-1))
    do isite=0,nsites-1
       iclass=lattice%myclass(isite,isite)
       do jclass=1,nlocclass
          if(iclass.eq.localtmp(jclass))exit
       enddo
       if(jclass>nlocclass)then
          nlocclass=nlocclass+1
          localtmp(jclass)=iclass
          Utmp(jclass)=hamilt%Uv(isite,isite)
          mutmp(jclass)=real(hamilt%hop(isite,isite))-hamilt%mu
       endif
       hamilt%mylocclass(isite)=jclass
    enddo
    hamilt%nlocclass=nlocclass
    allocate(hamilt%Uvalue(nlocclass),hamilt%muvalue(nlocclass))
    hamilt%Uvalue(1:nlocclass)=Utmp(1:nlocclass)
    hamilt%muvalue(1:nlocclass)=mutmp(1:nlocclass)
  end subroutine count_local_classes


end module DQMC_HAMILT

