module DQMC_GEOM_WRAP
  use DQMC_GEOM_PARAM
  use DQMC_HAMILT
  use DQMC_SYMM
  use DQMC_LATT
  use DQMC_RECLATT
  use DQMC_BONDS
  use DQMC_Cfg
  use DQMC_STRUCT
  implicit none


  type GeomWrap
     type(lattice_t)       :: Lattice
     type(recip_lattice_t) :: RecipLattice
     type(hamiltonian_t)   :: Hamilt
     type(bonds_t)         :: Bonds
     type(symm_operations) :: SymmOp
     type(pairing)         :: Pairs
   
     complex*16, pointer   :: FourierC(:,:)
  end type

  contains

  subroutine DQMC_Geom_Fill(gwrap,gfile)

    type(GeomWrap),    intent(inout)  :: gwrap
    character(len=30), intent(in)     :: gfile
    logical                           :: found, connected

    inquire(file=gfile,exist=found)
    if(found)then
     inquire(file=gfile,opened=connected)
     if(connected)then
      inquire(file=gfile,number=inpunit)
      rewind(inpunit)
     else
      inpunit=11
      open(file=gfile,unit=inpunit)
     endif
    else 
     call DQMC_Error("cannot open geom def file "//gfile, 0)
    endif

    call analyze_input
  
    !Initialize basic info about real space cluster
    call init_lattice(gwrap%Lattice)
    !Construct full lattice 
    call construct_lattice(gwrap%Lattice)
 
    !Initialize basic info about reciprocal lattice
    call init_recip_latt(gwrap%Lattice,gwrap%RecipLattice)
    !construct full lattice
    call construct_recip_lattice(gwrap%RecipLattice)

    !Construct Hamiltonian
    call construct_hamilt(gwrap%Hamilt,gwrap%Lattice,gwrap%RecipLattice)
    
    !Read point-symmetry (optional)
    call read_symm(gwrap%SymmOp)
    !pair-to-pair map of action of each symmetry in real space (also translations)
    call map_symm_lattice(gwrap%SymmOp,gwrap%Lattice)
    !point-to-point map of action of each symmetry in reciprocal space
    call map_symm_recip_lattice(gwrap%SymmOp,gwrap%RecipLattice)

    !Group pairs of lattice points into classes
    call construct_lattice_classes(gwrap%SymmOp,gwrap%Lattice)
    !Group k-points into classes
    call construct_recip_lattice_classes(gwrap%SymmOp,gwrap%RecipLattice)
    !Group hopping part of Hamiltonian in classes
    call count_hop_class(gwrap%Lattice,gwrap%Hamilt)
    !Group separately local classes for U and mu
    call count_local_classes(gwrap%Lattice,gwrap%Hamilt)

    !Read Bonds (optional input)   
    call read_bonds(gwrap%Bonds)
    if(gwrap%Bonds%initialized)then
     !bond-to-bond map of action of symmetry on bonds
     call map_symm_bonds(gwrap%Bonds,gwrap%SymmOp,gwrap%Lattice)
     !Group pairs of bonds into classes
     call construct_bond_classes(gwrap%Bonds,gwrap%SymmOp)
     !Map bonds throughout the entire lattice
     call construct_pairs(gwrap%Bonds,gwrap%Pairs,gwrap%Lattice)
     !pair-to-pair map of action of symmetry on bonds
     call map_symm_pairs(gwrap%Pairs, gwrap%SymmOp)
     !Group pairs of pairs into classes
     call construct_pair_classes(gwrap%Bonds,gwrap%Pairs)
    endif

    !Assign phase to each atom (optional input)
    call assign_phase(gwrap%Lattice)

    !Randomly dilute the lattice (optional input)
    call dilution(gwrap%Lattice,gwrap%Hamilt,gwrap%Pairs)

    !Write some info
    !call write_files(gwrap%Lattice,gwrap%RecipLattice,gwrap%Hamilt)

  end subroutine DQMC_Geom_Fill

 !---------------------------------------------------------------------!

  subroutine DQMC_Geom_Init(gwrap, S, cfg)    

    type(GeomWrap), intent(in)    :: gwrap   
    type(Struct),   intent(inout) :: S       ! Struct
    type(config),   intent(inout) :: cfg

    ! ... local scalar ...
    integer  :: n                ! Order of matrix T and D 
    integer  :: i, j             ! Loop iterator
    integer  :: idx, ineig, ic, ib              ! 
    real(wp), pointer  :: clab(:,:), tvalue(:) 
    integer, pointer   :: tmp(:,:)

    ! ... Executable ...
    S%checklist=.false.
  
    open(file='table.def',unit=25,status='unknown')

    n   = gwrap%Lattice%nsites
    S%nSite = n
    S%n_t=gwrap%Hamilt%nhopclass
    S%nClass=gwrap%Lattice%nclass
    S%nGroup=gwrap%Hamilt%nlocclass
    S%Name='General Geometry - Free Format'
    allocate(tmp(n,n))

    !Fill T
    tmp=0
    do i=0,n-1
     do ineig=1,gwrap%Hamilt%tnneig(i)
      j=gwrap%Hamilt%tneig(i,ineig)
      ic=gwrap%Hamilt%myhopclass(i,ineig)
      tmp(i+1,j+1)=ic
     enddo
    enddo
    call DQMC_CCS_Compress(n,-1, tmp, S%T)
    S%checklist(STRUCT_ADJ)=.true.
    
    !Fill B
    if(Found_Field(BONDS_F))then
     tmp=0
     S%n_b=gwrap%Pairs%nbond
     do ic=0,size(gwrap%Pairs%nbondv)-1
       do ib=1,gwrap%Pairs%nbondv(ic)
         i=gwrap%Pairs%bond_origin(ib,ic)
         j=gwrap%Pairs%bond_end(ib,ic)
         tmp(i+1,j+1)=gwrap%Pairs%bond_number(ib,ic)
      enddo
     enddo
     call DQMC_CCS_Compress(n,-1, tmp, S%B)

     !Store symmetry in S
     allocate(S%BC(S%n_b,S%n_b),S%BCsize(gwrap%Pairs%nclass_p))
     S%BC=gwrap%Pairs%myclass_p
     S%BCsize=gwrap%Pairs%class_size_p
     S%nBC=gwrap%Pairs%nclass_p
     S%checklist(STRUCT_BOND)=.true.

     !Waves
     S%nWave=gwrap%Pairs%nWave
     allocate(S%wlabel(S%nWave))
     S%wlabel(:)=gwrap%Pairs%wave_label(:)
     if(Found_Field(PAIRS_F))then
       call DQMC_reshape(S%n_b,S%nWave,S%W)
       do i=1,S%nWave
         S%W(1:S%n_b,i)=gwrap%Pairs%bond_wgt(i,1:S%n_b)
       enddo
       S%checklist(STRUCT_WAVE)=.true.
     endif
    endif

    deallocate(tmp)
    
    call DQMC_Reshape(n, n, S%D)
    call DQMC_Reshape(S%nClass, S%F)
    allocate(S%clabel(S%nClass))
    S%D(1:n,1:n)=gwrap%Lattice%myclass(0:n-1,0:n-1)
    S%F(:)=gwrap%Lattice%class_size(:)
    clab=>gwrap%Lattice%class_label
    do ic=1,S%nClass
     write(S%clabel(ic),'(i3,3(f8.4))') int(clab(ic,4)),(clab(ic,j),j=1,3)
    enddo
    S%checklist(STRUCT_CLASS)=.true.    

    call DQMC_Reshape(n, S%map)
    S%map(1:n)=gwrap%Hamilt%mylocclass(0:n-1)
    
    if(Found_field(PHASE_F))then 
     call DQMC_Reshape(n, S%P)
     S%P(1:n)=gwrap%Lattice%phase(0:n-1)
     S%checklist(STRUCT_PHASE)=.true.
    endif

    !Set variables that are otherwise read from main input
    call CFG_Set(cfg,"n",n)
    allocate(tvalue(S%n_t))
    tvalue(:)=dble(gwrap%Hamilt%tvalue(:))
    call CFG_Set(cfg,"t",S%n_t,tvalue)
    deallocate(tvalue)
    call CFG_Set(cfg,"U",S%nGroup,gwrap%Hamilt%Uvalue)
    call CFG_Set(cfg,"mu",S%nGroup,gwrap%Hamilt%muvalue)

    write(25,'(A,A)')'name = ',S%Name
    write(25,*)
    write(25,*)'nSite = ',n
    write(25,*)'dim ='
    write(25,*)
    write(25,*)'n_t =',gwrap%Hamilt%nhopclass
    write(25,*)'T = ',S%T%nnz
    call dqmc_CCS_Print(S%T,25)
    write(25,*)
    write(25,*)'nClass =',S%nClass
    write(25,*)'D =',n*n
    do i=1,n
      do j=1,n
        write(25,*)i,j,S%D(i,j)
      enddo
    enddo
    write(25,*)
    write(25,*)'clabel =',S%nClass
    do ic=1,S%nClass
      write(25,'(A)')S%clabel(ic)
    enddo
    if(Found_field(PHASE_F))then 
      write(25,*)
      write(25,*)'P =',n
      write(25,'(f10.5)')S%P(1:n)
    endif
    

    S%checklist(STRUCT_INIT)=.true.

  end subroutine DQMC_Geom_Init

 !---------------------------------------------------------------------!

  subroutine DQMC_Geom_Init_AntiPBC(gwrap, S, cfg)    

    type(GeomWrap), intent(in)    :: gwrap   
    type(Struct),   intent(inout) :: S       ! Struct
    type(config),   intent(inout) :: cfg

    ! ... local scalar ...
    integer  :: n                ! Order of matrix T and D 
    integer  :: i, j             ! Loop iterator
    integer  :: idx, ineig, ic, ib              ! 
    real(wp), pointer  :: clab(:,:), tvalue(:) 
    integer, pointer   :: tmp(:,:)

    ! ... Executable ...
    S%checklist=.false.
  
    open(file='table.def',unit=25,status='unknown')

    n   = gwrap%Lattice%nsites
    S%nSite = n
    S%n_t=gwrap%Hamilt%nhopclass
    S%nClass=gwrap%Lattice%nclass
    S%nGroup=gwrap%Hamilt%nlocclass
    S%Name='General Geometry - Free Format'
    allocate(tmp(n,n))

    !Fill T
    tmp=0
    do i=0,n-1
     do ineig=1,gwrap%Hamilt%tnneig(i)
      j=gwrap%Hamilt%tneig(i,ineig)
      ic=gwrap%Hamilt%myhopclass(i,ineig)
      tmp(i+1,j+1)=ic
     enddo
    enddo
    !CUSTOM MODIFICATION ASSUMING SQUARE LATTICE!!!!!
    !TO HANDLE ANTI-PBC!!!!!
    S%n_t=2
    ic=sqrt(dble(n))
    do i=ic,n,ic
      tmp(i-ic+1,i) = 2
      tmp(i,i-ic+1) = 2
    enddo
    j=n-ic
    do i=1,ic
      tmp(i,j+i)=2
      tmp(j+i,i)=2
    enddo
    call DQMC_CCS_Compress(n,-1, tmp, S%T)
    S%checklist(STRUCT_ADJ)=.true.
    
    !Fill B
    if(Found_Field(BONDS_F))then
     tmp=0
     S%n_b=gwrap%Pairs%nbond
     do ic=0,size(gwrap%Pairs%nbondv)-1
       do ib=1,gwrap%Pairs%nbondv(ic)
         i=gwrap%Pairs%bond_origin(ib,ic)
         j=gwrap%Pairs%bond_end(ib,ic)
         tmp(i+1,j+1)=gwrap%Pairs%bond_number(ib,ic)
      enddo
     enddo
     call DQMC_CCS_Compress(n,-1, tmp, S%B)

     !Store symmetry in S
     allocate(S%BC(S%n_b,S%n_b),S%BCsize(gwrap%Pairs%nclass_p))
     S%BC=gwrap%Pairs%myclass_p
     S%BCsize=gwrap%Pairs%class_size_p
     S%nBC=gwrap%Pairs%nclass_p
     S%checklist(STRUCT_BOND)=.true.

     !Waves
     S%nWave=gwrap%Pairs%nWave
     allocate(S%wlabel(S%nWave))
     S%wlabel(:)=gwrap%Pairs%wave_label(:)
     if(Found_Field(PAIRS_F))then
       call DQMC_reshape(S%n_b,S%nWave,S%W)
       do i=1,S%nWave
         S%W(1:S%n_b,i)=gwrap%Pairs%bond_wgt(i,1:S%n_b)
       enddo
       S%checklist(STRUCT_WAVE)=.true.
     endif
    endif

    deallocate(tmp)
    
    call DQMC_Reshape(n, n, S%D)
    call DQMC_Reshape(S%nClass, S%F)
    allocate(S%clabel(S%nClass))
    S%D(1:n,1:n)=gwrap%Lattice%myclass(0:n-1,0:n-1)
    S%F(:)=gwrap%Lattice%class_size(:)
    clab=>gwrap%Lattice%class_label
    do ic=1,S%nClass
     write(S%clabel(ic),'(i3,3(f8.4))') int(clab(ic,4)),(clab(ic,j),j=1,3)
    enddo
    S%checklist(STRUCT_CLASS)=.true.    

    call DQMC_Reshape(n, S%map)
    S%map(1:n)=gwrap%Hamilt%mylocclass(0:n-1)
    
    if(Found_field(PHASE_F))then 
     call DQMC_Reshape(n, S%P)
     S%P(1:n)=gwrap%Lattice%phase(0:n-1)
     S%checklist(STRUCT_PHASE)=.true.
    endif

    !Set variables that are otherwise read from main input
    call CFG_Set(cfg,"n",n)
    allocate(tvalue(S%n_t))
    tvalue(1)=dble(gwrap%Hamilt%tvalue(1))
    tvalue(2)=-tvalue(1)
    call CFG_Set(cfg,"t",S%n_t,tvalue)
    deallocate(tvalue)
    call CFG_Set(cfg,"U",S%nGroup,gwrap%Hamilt%Uvalue)
    call CFG_Set(cfg,"mu",S%nGroup,gwrap%Hamilt%muvalue)

    write(25,'(A,A)')'name = ',S%Name
    write(25,*)
    write(25,*)'nSite = ',n
    write(25,*)'dim ='
    write(25,*)
    write(25,*)'n_t =',gwrap%Hamilt%nhopclass
    write(25,*)'T = ',S%T%nnz
    call dqmc_CCS_Print(S%T,25)
    write(25,*)
    write(25,*)'nClass =',S%nClass
    write(25,*)'D =',n*n
    do i=1,n
      do j=1,n
        write(25,*)i,j,S%D(i,j)
      enddo
    enddo
    write(25,*)
    write(25,*)'clabel =',S%nClass
    do ic=1,S%nClass
      write(25,'(A)')S%clabel(ic)
    enddo
    if(Found_field(PHASE_F))then 
      write(25,*)
      write(25,*)'P =',n
      write(25,'(f10.5)')S%P(1:n)
    endif
    

    S%checklist(STRUCT_INIT)=.true.

  end subroutine DQMC_Geom_Init_AntiPBC


  !------------------------------------------------------------------------------------
  ! Remove some of the sites and update all relevant quantities. This routine
  ! assume that the number of classes is unchanged by this operation. Unless the
  ! dilution is really high this should always be the case.
  !------------------------------------------------------------------------------------
  subroutine dilution(lattice,hamilt,Pairs)
   type(lattice_t)   :: lattice
   type(Hamiltonian_t) :: hamilt
   type(pairing) :: Pairs
   integer :: nkept,isite,iold,jsite,jold,i,j,iclass,nsites2, &
              ios,iatom,icell,ibond,jbond,nsites,natom,ncell
   integer,allocatable :: kept_sites(:),old_to_new(:),tmp_bond_origin(:,:),&
                          tmp_bond_end(:,:),npairbondv(:)
   integer,pointer :: myclasstmp(:,:)
   complex*16, pointer :: hoptmp(:,:)
   real*8, pointer :: postmp(:,:),cartpostmp(:,:),Uvtmp(:,:),Jvtmp(:,:)
   real*8 rand,dilution_factor
   logical, allocatable :: to_be_kept(:)
   
   if(move_to_record(INPUT_FIELDS(DILUT_F),inpunit))then
    nsites=lattice%nsites
    natom=lattice%natom
    ncell=lattice%ncell
    allocate(kept_sites(nsites),old_to_new(0:nsites-1),to_be_kept(0:nsites-1))
    call random_seed
    nsites2=0
    do i=0,nsites-1
       call random_number(rand)
       if(rand<abs(dilution_factor))then
         to_be_kept(i)=.true.
         nsites2=nsites2+1
         kept_sites(nsites2)=i
       endif
    enddo
    if(nsites2==nsites)return  !nothing to do

    !nsites2=0
    !do 
    ! read(inpunit,*,iostat=ios)iatom,dilution_factor
    ! if(ios/=0)exit
    ! if(iatom>natom-1)stop 'Dilution applied to non-existing atom'
    ! do i=iatom,nsites-1,natom
    !   call random_number(rand)
    !   if(rand>dilution_factor)then
    !    nsites2=nsites2+1
    !    to_be_kept(i)=.true.
    !    kept_sites(nsites2)=i
    !   endif
    ! enddo
    !enddo
    !call isort(kept_sites,nsites2)

    write(*,*)'Sites to be removed'
    do isite=0,nsites-1
     if(.not.to_be_kept(isite))write(*,*)isite
    enddo
    allocate(myclasstmp(0:nsites2-1,0:nsites2-1),hoptmp(0:nsites2-1,0:nsites2-1),postmp(3,0:nsites2-1), &
     & cartpostmp(3,0:nsites2-1),Uvtmp(0:nsites2-1,0:nsites2-1),Jvtmp(0:nsites2-1,0:nsites2-1))
    !given old site label returns new. If old is eliminated returns -1.
    old_to_new(:)=-1
    do isite=0,nsites2-1
     old_to_new(kept_sites(isite+1))=isite
    enddo
    do isite=0,nsites2-1
     iold=kept_sites(isite+1)
     postmp(1:3,isite)=lattice%pos(1:3,iold)
     cartpostmp(1:3,isite)=lattice%cartpos(1:3,iold)
     write(*,'(i3,3f10.6)')isite,lattice%cartpos(1:3,iold)
     do jsite=0,nsites2-1
      jold=kept_sites(jsite+1)
      myclasstmp(isite,jsite)=lattice%myclass(iold,jold)
      hoptmp(isite,jsite)=hamilt%hop(iold,jold)
      Uvtmp(isite,jsite)=hamilt%Uv(iold,jold)
      Jvtmp(isite,jsite)=hamilt%Jv(iold,jold)
     enddo
    enddo
    !Dilute bonds for pairing
    if(Pairs%nbond>0)then
     allocate(tmp_bond_origin(Pairs%nbond,0:ncell-1),tmp_bond_end(Pairs%nbond,0:ncell-1),npairbondv(0:ncell-1))
     npairbondv(:)=0
     do icell=0,ncell-1
      do ibond=1,Pairs%nbond
       isite=Pairs%bond_origin(ibond,icell)
       jsite=Pairs%bond_end(ibond,icell)
       if(to_be_kept(isite).and.to_be_kept(jsite))then
        jbond=npairbondv(icell)+1
        npairbondv(icell)=jbond
        Pairs%bond_number(jbond,icell)=ibond
        tmp_bond_origin(jbond,icell)=old_to_new(isite)
        tmp_bond_end(jbond,icell)=old_to_new(jsite)
       endif
      enddo
     enddo
     Pairs%nbondv(:)=npairbondv(:)
     do icell=0,ncell-1
      jbond=npairbondv(icell)
      Pairs%bond_origin(1:jbond,icell)=tmp_bond_origin(1:jbond,icell)
      Pairs%bond_origin(jbond+1:Pairs%nbond,icell)=-1
      Pairs%bond_end(1:jbond,icell)=tmp_bond_end(1:jbond,icell)
      Pairs%bond_end(jbond+1:Pairs%nbond,icell)=-1
     enddo
     deallocate(tmp_bond_origin,tmp_bond_end,npairbondv)
    endif
    deallocate(lattice%myclass,hamilt%hop,lattice%pos,lattice%cartpos,hamilt%Uv,hamilt%Jv)
    deallocate(kept_sites,old_to_new,to_be_kept)
    !Update variable after site dilution
    lattice%myclass=>myclasstmp
    hamilt%hop=>hoptmp
    lattice%pos=>postmp
    lattice%cartpos=>cartpostmp
    hamilt%Uv=>Uvtmp
    hamilt%Jv=>Jvtmp
    lattice%nsites=nsites2
    lattice%class_size(:)=0
    do i=0,nsites2-1
     do j=0,nsites2-1
      iclass=lattice%myclass(i,j)
      lattice%class_size(iclass)=lattice%class_size(iclass)+1
     enddo
    enddo
    !Redefine neighbors and their classes
    call find_neighbors(hamilt)
    deallocate(hamilt%myhopclass,hamilt%tvalue)
    call count_hop_class(lattice,hamilt)
    deallocate(hamilt%mylocclass,hamilt%Uvalue,hamilt%muvalue)
    call count_local_classes(lattice,hamilt)
   endif
   
   contains
   
   subroutine isort(ivec,n)
   integer n,ivec(n),i,j,a,k
   do i=2,n
    a=ivec(i)
    do j=1,i-1
     if(a<ivec(j))exit
    enddo
    do k=i-1,j,-1
     ivec(k+1)=ivec(k)
    enddo
    ivec(j)=a
   enddo
  end subroutine
  
  end subroutine dilution



  !---------------------------------------------------------------
  ! Fill the matrix of Fourier coefficients
  !--------------------------------------------------------------
  subroutine DQMC_Fill_FourierC(gwrap)
  type(GeomWrap),target :: gwrap
  integer               :: nt,nk,i,ii,j
  real*8, pointer       :: tr(:,:),kpts(:,:)
  integer, pointer      :: indx(:)

  !initialize
  tr     => gwrap%lattice%translation
  kpts   => gwrap%RecipLattice%klist
  indx   => gwrap%RecipLattice%class_repr_k
  nt = gwrap%lattice%ncell
  nk = gwrap%RecipLattice%nclass_k
  allocate(gwrap%FourierC(nt,nk))

  !compute
  do i= 0, nt-1
   ii=i+1
   do j= 1, nk
    gwrap%FourierC(ii,j)=exp(im*sum(tr(:,i)*kpts(indx(j),:)))
   enddo
  enddo

  end subroutine DQMC_Fill_FourierC


!  !-----------------------------------------------------------------------
!  ! Write extended files 
!  !-----------------------------------------------------------------------
!  subroutine write_files(lattice,recip_lattice,hamilt)
!  type(lattice_t)   :: lattice
!  type(Hamiltonian_t) :: hamilt
!  type(recip_lattice_t) :: recip_lattice
!  integer isite,ineig,ihc,jsite,nsites,nhopclass,nlocclass
!  character*50 inpfname,fname
!  
!  nsites=lattice%nsites
!  nhopclass=hamilt%nhopclass
!  nlocclass=hamilt%nlocclass
!  
!  inquire(unit=inpunit,name=inpfname)
!  inpfname=adjustl(inpfname)
!  
!  fname=trim(inpfname)//'.neigh2'
!  open(unit=25,file=fname)
!  write(25,*) nsites,nhopclass
!  do isite=0,nsites-1
!   do ineig=1,hamilt%tnneig(isite)
!    write(25,*) isite+1,hamilt%tneig(isite,ineig)+1,hamilt%myhopclass(isite,ineig)
!   enddo
!  enddo
!  close(25)
!  
!  fname=trim(inpfname)//'.tmap2'
!  open(unit=25,file=fname)
!  write(25,*) nhopclass
!  do ihc=1,nhopclass
!   write(25,*)ihc,hamilt%tvalue(ihc)
!  enddo
!  close(25)
!  
!  fname=trim(inpfname)//'.local2'
!  open(unit=25,file=fname)
!  write(25,*) nsites,nlocclass
!  do isite=0,nsites-1
!   write(25,*) isite+1,hamilt%mylocclass(isite)
!  enddo
!  close(25)
!  
!  fname=trim(inpfname)//'.Umap2'
!  open(unit=25,file=fname)
!  write(25,*) nlocclass
!  do ihc=1,nlocclass
!   write(25,*)ihc,hamilt%Uvalue(ihc)
!  enddo
!  close(25)
!  
!  fname=trim(inpfname)//'.mumap2'
!  open(unit=25,file=fname)
!  write(25,*) nlocclass
!  do ihc=1,nlocclass
!   write(25,*)ihc,hamilt%muvalue(ihc)
!  enddo
!  close(25)
!  
!  fname=trim(inpfname)//'.class2'
!  open(unit=25,file=fname)
!  write(25,*) nsites,lattice%nclass
!  do isite=0,nsites-1
!   do jsite=0,nsites-1
!    write(25,*) isite+1,jsite+1,lattice%myclass(isite,jsite)
!   enddo
!  enddo
!  close(25)
!  
!  fname=trim(inpfname)//'.phase2'
!  open(unit=25,file=fname)
!  write(25,*) nsites
!  do isite=0,nsites-1
!   write(25,*)isite+1,lattice%susc_phase(isite)
!  enddo
!  close(25)
!  
!  fname=trim(inpfname)//'.kspace'
!  open(unit=25,file=fname)
!  write(25,*) recip_lattice%nclass_k,recip_lattice%nkpts
!  do ihc=1,recip_lattice%nkpts
!   write(25,'(i3,3(1x,f10.5),i3)')ihc,(recip_lattice%klist(ihc,isite),isite=1,3),recip_lattice%myclass_k(ihc)
!  enddo
!  close(25)
!  
!  end subroutine


end module DQMC_GEOM_WRAP
