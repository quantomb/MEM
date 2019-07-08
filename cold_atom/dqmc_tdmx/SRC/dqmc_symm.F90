module DQMC_SYMM

  use DQMC_GEOM_PARAM

  implicit none

  type :: symm_operations
     integer               :: nsymm                 !number of symmetris compatible with supercell
     integer               :: ntransl               !number of translations
     integer               :: ntotsymm              !number of symmetry read from input
     integer, pointer      :: map_symm(:,:)         !action of symmetry on site (0:nsites-1, nsymm)
     integer, pointer      :: map_symm_k(:,:)       !action of symmetry on k-point (nkpts, nsymm)
     integer, pointer      :: map_symm_b(:,:)       !action of symmetry on bond (ntotbond, nsymm)
     integer, pointer      :: map_symm_p(:,:)       !action of symmetry on pair (nbond, nsymm)
     integer, pointer      :: translback(:)         !number of translation mapping a site into its untranslated image (0:nsites-1) 
     integer, pointer      :: translate(:,:)        !action of translations on site (0:nsites-1,0:ntransl-1)
     integer, pointer      :: valid_symm(:)         !label of valid symmetry operation(nsymm)

     real*8, pointer       :: symmangle(:)          !Angle of rotation of symmetry operation (ntotsymm)
     real*8, pointer       :: symmpoint(:,:)        !point which symmetry operation goes through (rdim,ntotsymm)
     real*8, pointer       :: symmaxis(:,:)         !direction of the axis (rdim,ntotsymm)

     character*1, pointer  :: symmlabel(:)          !label of operation: C, D, I (ntotsymm)

     logical               :: initialized
     logical               :: lattice_mapped
     logical               :: recip_lattice_mapped
     logical               :: bonds_mapped
     logical               :: TimeRev

  end type symm_operations

contains



  !------------------------------------------------------------------------------
  ! Read the point-symmetries from input.
  ! Rotation-axis : CN x y z x1 y1 z1
  !  2\pi/N is the rotation angle
  !  x y z specify a point belonging to the axis in cartesian coordinates
  !  x1 y1 z1 specify the axis direction in cartesian coordinates
  ! Mirror plane : D x y z x1 y1 z1
  !  x y z specify a point belonging to the plane in cartesian coordinate
  !  x1 y1 z1 specify the direction normal to the plane in cartesian coordinates
  ! Inversion : I x y z
  !  x y z specify position of inversion point in cartesian coordinates
  !------------------------------------------------------------------
  subroutine read_symm(SymmOp)
    integer isymm,axis_order,ios1,ios2,i,nsymm
    real*8 xpoint(3),xaxis(3),xnorm
    character*50 string
    character*1 label
    logical dum
    character*1, pointer :: symmlabel(:)
    real*8, pointer      :: symmangle(:),symmpoint(:,:),symmaxis(:,:)
    type(symm_operations),intent(out) :: SymmOp
    nsymm=count_symmetry()
    SymmOp%initialized=.false.
    SymmOp%nsymm=nsymm
    SymmOp%ntotsymm=nsymm
    SymmOp%TimeRev=.true.
    if(nsymm/=0)then
       allocate(symmangle(nsymm),symmpoint(3,nsymm),symmaxis(3,nsymm),symmlabel(nsymm))
       dum=move_to_record(INPUT_FIELDS(SYMM_F),inpunit)
       do isymm=1,nsymm
          read(inpunit,'(A)')string
          read(string,'(A1)')label
          if(label=='C'.or.label=='c')then
             read(string,'(A1,i1)',iostat=ios1)label,axis_order
             read(string(3:50),*,iostat=ios2)(xpoint(i),i=1,3),(xaxis(i),i=1,3)
             if(ios1.ne.0.or.ios2.ne.0)then
                write(*,*)'Problem reading axis symmetry. Stop. Line:',isymm
                stop
             endif
             symmlabel(isymm)=label(1:1)
             symmangle(isymm)=4.d0*acos(0.d0)/axis_order
             symmpoint(:,isymm)=xpoint(:)
             xnorm=sqrt(sum(xaxis(:)**2))
             symmaxis(:,isymm)=xaxis(:)/xnorm
          elseif(label=='D'.or.label=='d')then
             read(string(3:50),*,iostat=ios2)(xpoint(i),i=1,3),(xaxis(i),i=1,3)
             if(ios2.ne.0)then
                write(*,*)'Problem reading plane symmetry. Stop. Line:',isymm
                stop
             endif
             symmlabel(isymm)=label(1:1)
             symmangle(isymm)=0.d0
             symmpoint(:,isymm)=xpoint(:)
             xnorm=sqrt(sum(xaxis(:)**2))
             symmaxis(:,isymm)=xaxis(:)/xnorm
          elseif(label=='I'.or.label=='i')then
             read(string(3:50),*,iostat=ios2)(xpoint(i),i=1,3)
             if(ios2.ne.0)exit
             symmlabel(isymm)=label(1:1)
             symmangle(isymm)=0.d0
             symmpoint(:,isymm)=xpoint(:)
             symmaxis(:,isymm)=0.d0
             SymmOp%TimeRev=.false.
          else
             write(*,*)'Unknown label for symmetry. Stop'
             stop
          endif
       enddo
       !Fill SymmOp
       SymmOp%initialized=.true.
       SymmOp%nsymm=nsymm
       SymmOp%ntotsymm=nsymm
       SymmOp%symmlabel=>symmlabel
       SymmOp%symmaxis=>symmaxis
       SymmOp%symmpoint=>symmpoint
       SymmOp%symmangle=>symmangle
    endif
  end subroutine read_symm



  !--------------------------------------------------------
  !Count how many point-symmetries are specified in input
  !--------------------------------------------------------
  !integer function count_symmetry
  function count_symmetry()
    integer :: msymm, count_symmetry
    character*1 ::label
    logical :: ldum
    msymm=0
    if(Found_Field(SYMM_F))then
       ldum=move_to_record(INPUT_FIELDS(SYMM_F),inpunit)
       do
          read(inpunit,'(A1)')label
          if(   label=='C'.or.label=='c' &
               & .or.label=='D'.or.label=='d' &
               & .or.label=='I'.or.label=='i' )then
             msymm=msymm+1
          else
             exit
          endif
       enddo
    else
       write(*,*)'No symmetries have been specified'
    endif
    count_symmetry=msymm
  end function count_symmetry




  !----------------------------------------------------------------------
  ! Apply point-symmetry operation. label specifis the kind of
  ! symmetry (C, D or I). Point and axis specify the "position" of
  ! the symmetry operation (see comments to read_symmetry). set is a set
  ! of 3D vectors in cartesian coordinates upon which the symmetry acts.
  ! newset returns the transformed set.
  !-----------------------------------------------------------------------
  subroutine apply_point_symm(label,point,axis,theta,set,new_set,nset,reciprocal)
    integer,intent(in) :: nset
    real*8,intent(in) :: point(3),axis(3),theta,set(3*nset)
    logical, intent(in) :: reciprocal
    real*8,intent(out) :: new_set(3,nset)
    integer irefl,iset
    real*8 trans_set(3,nset),d,a1,a2,rot1(3,3),rot2(3,3),rot3(3,3),globrot(3,3), &
         &cost,sint,dpoint(3),dset(3,nset),dset2(nset,3)
    integer:: i,j,alp,bet,gam,del
    character*1 label
    if(reciprocal)then
       dpoint(:)=0.d0
       dset=reshape(set,shape=(/3,nset/),order=(/2,1/))
    else
       dpoint(:)=point(:)
       dset=reshape(set,shape=(/3,nset/),order=(/1,2/))
    endif
    if(label=='i'.or.label=='I')then
       do iset=1,nset
          new_set(:,iset)=2*dpoint(:)-dset(:,iset)
       enddo
    else
       !Define matrix that line up z-axis with rotation axis
       if(label=='c'.or.label=='C')then 
          irefl=1
       else
          irefl=-1
       endif
       cost=cos(theta); sint=sin(theta)
       do iset=1,nset
          trans_set(:,iset)=dset(:,iset)-dpoint(:)
       enddo
       d=sqrt(axis(1)**2+axis(2)**2)
       if(d>1.d-6)then
          a1=axis(1)/d; a2=axis(2)/d
       else
          a1=1.d0; a2=0.d0
       endif
       !First rotate reference frame around z : R1
       rot1(1,1)= a1       ; rot1(1,2)= a2       ; rot1(1,3)=0.d0
       rot1(2,1)=-a2       ; rot1(2,2)= a1       ; rot1(2,3)=0.d0
       rot1(3,1)= 0.d0     ; rot1(3,2)= 0.d0     ; rot1(3,3)=1.d0
       !then Rotate reference frame around y : R2
       rot2(1,1)= axis(3)  ; rot2(1,2)= 0.d0     ; rot2(1,3)= -d
       rot2(2,1)= 0.d0     ; rot2(2,2)= 1.d0     ; rot2(2,3)= 0.d0
       rot2(3,1)= d        ; rot2(3,2)= 0.d0     ; rot2(3,3)= axis(3)
       !Finally rotate around z or reflect perpendicularly to xy (when irefl=-1): R3
       rot3(1,1)= cost      ; rot3(1,2)=-sint    ; rot3(1,3)=0.d0
       rot3(2,1)= sint      ; rot3(2,2)= cost    ; rot3(2,3)=0.d0
       rot3(3,1)= 0.d0      ; rot3(3,2)= 0.d0    ; rot3(3,3)=1.d0*irefl
       !compute globrot as R1^-1 R2^-1 R3 R2 R1 (undoing frame rotations)
       globrot(:,:)=0.d0
       do i  =1,3
       do j  =1,3
       do alp=1,3
       do bet=1,3
       do gam=1,3
       do del=1,3
             globrot(i,j)=globrot(i,j)+rot1(alp,i)*rot2(bet,alp)*rot3(bet,gam)*rot2(gam,del)*rot1(del,j)
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo
       !Rotate
       new_set(:,:)=0.d0
       do iset=1,nset
          do i=1,3
             do j=1,3
                new_set(i,iset)=new_set(i,iset)+globrot(i,j)*trans_set(j,iset)
             enddo
          enddo
       enddo
       !translate back to the origin
       do iset=1,nset
          new_set(:,iset)=new_set(:,iset)+dpoint(:)
       enddo
    endif
    if(reciprocal)then
       !Fancy way to get transpose
       dset2=reshape(new_set,shape=(/nset,3/),order=(/2,1/))
       new_set=reshape(dset2,shape=(/3,nset/))
    endif
  end subroutine apply_point_symm




  !---------------------------------------------------------------------
  ! Map the action of the point-operation into the vector map_symm:
  !map_symm(ifrom,isymm) returns the orbital on which orbital ifrom
  !is transformed by isymm. Symmetries that are incompatible with the
  !supercell are discarded. msymm is the number of compatible point
  !symmetries.
  ! Construct translate(ifrom,itrans). Analogous to map_symm but returns
  !the orbital to which ifrom is mapped by translation itrans. 
  ! Construct translback(ifrom). It returns the translation that applied
  !to ifrom reduces the latter to the equivalent orbital inside the
  !primitive cell.
  !---------------------------------------------------------------------
  subroutine map_symm_lattice(SymmOp,lattice)
    use DQMC_LATT
    type(symm_operations) :: SymmOp
    type(lattice_t),intent(in)   :: lattice
    integer :: i,msymm,iat,jat,ii,valid_symm(SymmOp%ntotsymm),it,itl,is,ntransl,natom, &
         & tmp_symm(0:lattice%nsites-1,SymmOp%ntotsymm),istart,ipat,jpat,nsites,ntotsymm,ndim
    real*8 newpos(rdim,0:lattice%nsites-1),diff(rdim),projsc(rdim),invscc(rdim,rdim)
    logical mapped(0:lattice%nsites-1,SymmOp%ntotsymm),mappedt(0:lattice%nsites-1,0:lattice%ncell-1)
    character*1, pointer :: symmlabel(:)
    real*8, pointer      :: symmangle(:),symmpoint(:,:),symmaxis(:,:)
    !Assign pointers
    symmlabel=>SymmOp%symmlabel; symmangle=>SymmOp%symmangle
    symmpoint=>SymmOp%symmpoint; symmaxis=>SymmOp%symmaxis
    !Initializa local variables
    ntransl=lattice%ncell
    nsites=lattice%nsites
    natom=lattice%natom
    ntotsymm=SymmOp%ntotsymm
    ndim=lattice%ndim
    call get_inverse(lattice%scc,invscc)
    !Apply symmetry operations to sites...
    mapped(:,:)=.false.
    do i=1,ntotsymm
       call apply_point_symm(symmlabel(i),symmpoint(:,i),symmaxis(:,i),symmangle(i),lattice%cartpos,newpos,nsites,.false.)
       !write(*,*)
       do iat=0,nsites-1
          ipat=mod(iat,natom)
          do jat=0,nsites-1
             jpat=mod(iat,natom)
             if(mapped(jat,i).or.lattice%olabel(jpat)/=lattice%olabel(ipat))cycle
             diff(:)=newpos(:,iat)-lattice%cartpos(:,jat)
             do ii=1,3
                projsc(ii)=sum(diff(:)*invscc(ii,:))
             enddo
             !only project inside the supercell along the "extended" dimensions
             projsc(1:ndim)=projsc(1:ndim)-nint(projsc(1:ndim))
             if(sum(projsc(1:rdim)**2)<toll)then
                !write(*,'(i3,i3,3f12.6)')iat,jat,diff(:)
                tmp_symm(iat,i)=jat
                mapped(jat,i)=.true.
                exit
             endif
          enddo
       enddo
    enddo
    !Check whether symmetry is compatible with supercell 
    msymm=0
    do i=1,ntotsymm
       do iat=0,nsites-1
          if(.not.mapped(iat,i))exit
       enddo
       if(iat==nsites)then
          msymm=msymm+1
          valid_symm(msymm)=i
       endif
    enddo
    if(msymm>0)allocate(SymmOp%map_symm(0:nsites-1,msymm),SymmOp%valid_symm(msymm))
    do i=1,msymm
       SymmOp%valid_symm(i)=valid_symm(i)
       SymmOp%map_symm(:,i)=tmp_symm(:,valid_symm(i))
    enddo
    !Deal with translation
    SymmOp%ntransl=ntransl
    if(ntransl>0)allocate(SymmOp%translate(0:nsites-1,0:ntransl-1),SymmOp%translback(0:nsites-1))
    mappedt(:,:)=.false.
    do it=0,ntransl-1
       itl=it*natom
       do iat=0,nsites-1
          newpos(:,iat)=lattice%cartpos(:,iat)+lattice%translation(:,it)
       enddo
       do iat=0,nsites-1
          istart=mod(iat,natom)
          do jat=istart,nsites-1,natom
             if(mappedt(jat,it))cycle
             diff(:)=newpos(:,iat)-lattice%cartpos(:,jat)
             do ii=1,3
                projsc(ii)=sum(diff(:)*invscc(ii,:))
             enddo
             projsc(1:ndim)=projsc(1:ndim)-nint(projsc(1:ndim))
             if(sum(projsc(1:rdim)**2)<toll)then
                SymmOp%translate(iat,it)=jat
                mappedt(jat,it)=.true.
                exit
             endif
          enddo
       enddo
    enddo
    do is=0,nsites-1
       do it=0,ntransl-1
          iat=SymmOp%translate(is,it)
          if(iat<natom)then
             SymmOp%translback(is)=it
             exit
          endif
       enddo
    enddo
    SymmOp%nsymm=msymm
    !Write out symmetry mappings
    write(25,*)'NSYMM',ntotsymm
    do i=1,ntotsymm
       write(25,*)'Simmetry',valid_symm(i),Symmlabel(valid_symm(i))
       do iat=0,nsites-1
          write(25,*)iat,'->',SymmOp%map_symm(iat,i)
       enddo
    enddo
    do it=0,ntransl-1
       write(25,*)'translation',it
       do iat=0,nsites-1
          write(25,*)iat,'->',SymmOp%translate(iat,it)
       enddo
    enddo
    write(25,*)'Translback'
    do iat=0,nsites-1
       write(25,*)iat,'->',SymmOp%translback(iat)
    enddo
    SymmOp%lattice_mapped=.true.
  end subroutine map_symm_lattice




  !--------------------------------------------------------------------------------
  !k-space symmetry
  !--------------------------------------------------------------------------------
  subroutine map_symm_recip_lattice(SymmOp,recip_lattice)
    use DQMC_RECLATT
    type(symm_operations)              :: SymmOp
    type(recip_lattice_t),intent(in)   :: recip_lattice
    integer                            :: i,j,ii,nsymm,ndim,ik,jk
    real*8                             :: diff(rdim),projsc(rdim),invkc(rdim,rdim), &
         newklist(recip_lattice%nkpts,rdim),zerovec(rdim)
    logical                            :: mapped(recip_lattice%nkpts)
    character*1, pointer               :: symmlabel(:)
    real*8, pointer                    :: symmangle(:),symmpoint(:,:),symmaxis(:,:)

    if(.not.SymmOp%lattice_mapped)stop'Need to map lattice symmetries before recip lattice ones'

    !Assign pointers
    symmlabel=>SymmOp%symmlabel; symmangle=>SymmOp%symmangle
    symmpoint=>SymmOp%symmpoint; symmaxis=>SymmOp%symmaxis
    nsymm=SymmOp%nsymm

    !Check whether we need to additionally apply time-reversal symmetry
    if(.not.SymmOp%TimeRev)then
       do i=1,nsymm
          j=SymmOp%valid_symm(i)
          if(SymmOp%symmlabel(j)=='i'.or.SymmOp%symmlabel(j)=='I')exit
       enddo
       if(i>nsymm)SymmOp%TimeRev=.true.
    endif
    if(SymmOp%TimeRev)nsymm=nsymm+1

    ndim=recip_lattice%ndim
    call get_inverse(recip_lattice%kc,invkc)
    if(nsymm>0)allocate(SymmOp%map_symm_k(recip_lattice%nkpts,nsymm))
    do j=1,SymmOp%nsymm
       i=SymmOp%valid_symm(j)
       mapped(:)=.false.
       call apply_point_symm(symmlabel(i),symmpoint(:,i),symmaxis(:,i),symmangle(i),&
            recip_lattice%klist,newklist,recip_lattice%nkpts,.true.)
       do ik=1,recip_lattice%nkpts
          do jk=1,recip_lattice%nkpts
             if(mapped(jk))cycle
             diff(:)=newklist(ik,:)-recip_lattice%klist(jk,:)
             do ii=1,rdim
                projsc(ii)=sum(diff(:)*invkc(:,ii))
             enddo
             !only project inside the supercell along the "extended" dimensions
             projsc(1:ndim)=projsc(1:ndim)-nint(projsc(1:ndim))
             if(sum(projsc(1:rdim)**2)<toll)then
                SymmOp%map_symm_k(ik,j)=jk
                mapped(jk)=.true.
                exit
             endif
          enddo
          if(jk>recip_lattice%nkpts)stop'Problem with symmetry in k-space'
       enddo
    enddo

    if(SymmOp%TimeRev)then
       mapped(:)=.false.
       zerovec(1:3)=(/0.d0, 0.d0, 0.d0/)
       call apply_point_symm('I', zerovec, zerovec, 0.d0,&
            recip_lattice%klist,newklist,recip_lattice%nkpts,.true.)
       do ik=1,recip_lattice%nkpts
          do jk=1,recip_lattice%nkpts
             if(mapped(jk))cycle
             diff(:)=newklist(ik,:)-recip_lattice%klist(jk,:)
             do ii=1,rdim
                projsc(ii)=sum(diff(:)*invkc(:,ii))
             enddo
             !only project inside the supercell along the "extended" dimensions
             projsc(1:ndim)=projsc(1:ndim)-nint(projsc(1:ndim))
             if(sum(projsc(1:rdim)**2)<toll)then
                SymmOp%map_symm_k(ik,nsymm)=jk
                mapped(jk)=.true.
                exit
             endif
          enddo
          if(jk>recip_lattice%nkpts)stop'Problem with symmetry in k-space'
       enddo
    endif
    SymmOp%recip_lattice_mapped=.true.
  end subroutine map_symm_recip_lattice




  !---------------------------------------------------------------------------------
  ! Given the list of bonds read from input, this routines does the following:
  ! 1) complete the list with the bonds which are equivalent, by symmetry,
  !    to those specified in input.
  ! 2) Create a mapping (map_symm_b) that, given, bond "b" and a symmetry operation "s",
  !    returns a the bond on which "b" is mapped by "s".
  !---------------------------------------------------------------------------------
  subroutine map_symm_bonds(Bonds,SymmOp,Lattice)
    use DQMC_BONDS
    type(symm_operations),intent(inout) :: SymmOp
    type(bonds_t),intent(inout) :: Bonds
    type(lattice_t),intent(in) :: Lattice
    integer :: nbclass,ibond,iat,jat,ntotpair,natom,nsites,bcl,newlabel,jbond,it, &
         & isymm,ndim
    integer, allocatable :: class(:),tag(:),pair_origin(:),pair_target(:),pair_label(:)
    integer, pointer :: bond_origin(:),bond_target(:),map_symm(:,:)
    real*8 :: invscc(rdim,rdim),proj(rdim)
    real*8, allocatable :: xxpair(:,:)
    logical, allocatable :: bond_on(:,:)

    if(.not.lattice%analyzed)stop'Need to analyze lattice before mapping bonds'
    if(Bonds%ntotbond==0)return

    natom=Lattice%natom; nsites=Lattice%nsites; ndim=Lattice%ndim
    call get_inverse(lattice%scc,invscc)
    nbclass=0
    newlabel=maxval(Bonds%bond_label)
    allocate(class(Bonds%ntotbond),tag(Bonds%ntotbond),bond_on(0:natom-1,0:nsites-1),bond_target(Bonds%ntotbond))
    bond_origin=>Bonds%bond_origin
    bond_on(:,:)=.false.
    do ibond=1,Bonds%ntotbond
       iat=bond_origin(ibond)
       jat=hoptowho(iat,Bonds%xxbond(:,ibond),Bonds%bond_target(ibond),Lattice)
       bond_target(ibond)=jat
       bond_on(iat,jat)=.true.
       class(ibond)=Lattice%myclass(iat,jat)
       !See if bond belongs to an already found class
       do jbond=1,ibond-1
          if(class(ibond)==class(jbond))exit
       enddo
       if(jbond==ibond)then
          !if not, create a new class
          nbclass=nbclass+1
          tag(nbclass)=class(ibond)
       endif
    enddo
    !Include all bonds which were left out but that are equivalent by symmetry...
    ntotpair=2*natom*nsites
    allocate(pair_label(ntotpair),pair_origin(ntotpair),pair_target(ntotpair),xxpair(3,ntotpair))
    ntotpair=Bonds%ntotbond
    pair_origin(1:ntotpair)=bond_origin(1:ntotpair)
    pair_target(1:ntotpair)=bond_target(1:ntotpair)
    pair_label(1:ntotpair)=Bonds%bond_label(1:ntotpair)
    xxpair(:,1:ntotpair)=Bonds%xxbond(:,1:ntotpair)
    !... But only if there was no PAIR field specified in input
    if(.not.Found_Field(PAIRS_F))then
       do iat=0,natom-1
          do jat=0,nsites-1
             if(bond_on(iat,jat))cycle
             do bcl=1,nbclass
                if(tag(bcl)==Lattice%myclass(iat,jat))then
                   newlabel=newlabel+1
                   ntotpair=ntotpair+1
                   pair_origin(ntotpair)=iat
                   pair_target(ntotpair)=jat
                   pair_label(ntotpair)=newlabel
                   xxpair(:,ntotpair)=lattice%cartpos(:,jat)-lattice%cartpos(:,iat)
                   !make xxpair as small as possible
                   do it=1,rdim
                      proj(it)=sum(xxpair(:,ntotpair)*invscc(it,:))
                   enddo
                   proj(1:ndim)=proj(1:ndim)-nint(proj(1:ndim))
                   do it=1,rdim
                      xxpair(it,ntotpair)=sum(lattice%scc(it,:)*proj(:))
                   enddo
                   bond_on(iat,jat)=.true.
                   if(iat/=jat)then
                      ntotpair=ntotpair+1
                      pair_label(ntotpair)=-newlabel
                      it=SymmOp%translback(jat) 
                      pair_origin(ntotpair)=SymmOp%translate(jat,it)
                      pair_target(ntotpair)=SymmOp%translate(iat,it)
                      xxpair(:,ntotpair)=-xxpair(:,ntotpair-1)
                      bond_on(SymmOp%translate(jat,it),SymmOp%translate(iat,it))=.true.
                   endif
                   exit
                endif
             enddo
          enddo
       enddo
    endif
    deallocate(class,tag,bond_target,bond_on)
    if(Bonds%ntotbond/=ntotpair)then
       !reload Bonds (new set completed with newly found bonds)
       deallocate(Bonds%bond_origin,Bonds%xxbond,Bonds%bond_label,Bonds%bond_target)
       allocate(Bonds%bond_origin(ntotpair),Bonds%bond_target(ntotpair), &
            Bonds%xxbond(3,ntotpair),Bonds%bond_label(ntotpair))
       Bonds%ntotbond=ntotpair
       Bonds%bond_origin(1:ntotpair)=pair_origin(1:ntotpair)
       Bonds%bond_target(1:ntotpair)=mod(pair_target(1:ntotpair),natom)
       Bonds%xxbond(:,1:ntotpair)=xxpair(:,1:ntotpair)
       Bonds%bond_label(1:ntotpair)=pair_label(1:ntotpair)
    endif
    allocate(map_symm(ntotpair,SymmOp%nsymm))
    !Store how a bond transforms under point symmetry
    do isymm=1,SymmOp%nsymm
       do ibond=1,ntotpair
          iat=SymmOp%map_symm(pair_origin(ibond),isymm)
          jat=SymmOp%map_symm(pair_target(ibond),isymm)
          it=SymmOp%translback(iat) 
          iat=SymmOp%translate(iat,it)
          jat=SymmOp%translate(jat,it)
          do jbond=1,ntotpair
             if(pair_origin(jbond)==iat.and.pair_target(jbond)==jat)exit
          enddo
          if(jbond>ntotpair)stop 'Symmetry analysis : Cannot find equivalent bond'
          map_symm(ibond,isymm)=jbond
       enddo
    enddo
    SymmOp%map_symm_b=>map_symm
    SymmOp%bonds_mapped=.true.
    write(*,*)
    write(*,*)'Bonds (Set completed using symmetry)'
    do ibond=1,Bonds%ntotbond
       write(*,'(3i4,3f12.7)')ibond,Bonds%bond_label(ibond),Bonds%bond_origin(ibond),Bonds%xxbond(1:rdim,ibond)
    enddo
    !write(*,*)
    !write(*,*)'BOND MAPPING'
    !do isymm=1,SymmOp%nsymm
    ! write(*,*)'Symmetry',isymm
    ! do ibond=1,Bonds%ntotbond
    !  write(*,*)ibond,'-->',map_symm(ibond,isymm)
    ! enddo
    !enddo
    write(*,*)'====================================================================='
    deallocate(pair_origin,pair_target,pair_label,xxpair)
    Bonds%analyzed=.true.
  end subroutine map_symm_bonds




  !--------------------------------------------------------------------------------
  ! Map symmetry for pairs.
  !--------------------------------------------------------------------------------
  subroutine map_symm_pairs(Pairs, SymmOp)

    use DQMC_BONDS
    type(symm_operations),intent(in) :: SymmOp
    type(pairing),intent(inout) :: Pairs
    integer :: isymm, ib, jb, newjb
    integer, pointer :: map_symm_p(:,:)

    map_symm_p => SymmOp%map_symm_p
    allocate(map_symm_p(Pairs%nbond,SymmOp%nsymm))
    !allocate(SymmOp%map_symm_p(Pairs%nbond,SymmOp%nsymm))
    do isymm=1,SymmOp%nsymm
       do ib=1,Pairs%nbond
          jb=Pairs%bond_map(ib)
          newjb=SymmOp%map_symm_b(jb,isymm)
          jb=Pairs%pair_map(newjb)
          if(jb==0)stop 'Symmetry analysis : Cannot find equivalent pair'
          !SymmOp%map_symm_p(ib,isymm)=jb
          map_symm_p(ib,isymm)=jb
       enddo
    enddo
  end subroutine map_symm_pairs




  !------------------------------------------------------------------------------
  !  Construct myclass(i,j). Given to sites i and j (not necessarily different)
  ! returns the class to which they belong. A class contains pairs of orbitals
  ! that, because of symmetry, are going to have identical pair-correlation
  ! functions.
  !  Returns nclass, the number of classes, and class_size(iclass), the 
  ! number of pairs inside class iclass.
  !  Returns class_label. This is a 4-components array. The first three are the
  ! cartesian separation of the two orbitals in the pair. The last component is
  ! the number of the atom inside the primitive cell that belongs to the pair.
  !------------------------------------------------------------------------------
  subroutine construct_lattice_classes(SymmOp,lattice)
    use DQMC_LATT
    integer               :: i,it,ip,is,j,isymm,iclass,istart,csize,csizenew,itransl,&
         &jclass,idj,id,mclass, ip_transl,is_transl,ip2,is2,jstart,&
         &nclass,nsites,natom,nsymm,ntransl
    integer,allocatable   ::  patom(:,:),satom(:,:),csizev(:)
    integer, pointer      :: myclass(:,:)
    type(symm_operations) :: SymmOp
    type(lattice_t)       :: lattice

    if(.not.SymmOp%lattice_mapped)stop'Need to map symmetries over lattice before classes'

    !initialize local variables
    natom=lattice%natom
    nsites=lattice%nsites
    nsymm=SymmOp%nsymm
    ntransl=SymmOp%ntransl
    !allocate internal arrays
    allocate(patom(natom*nsites,natom*nsites),satom(natom*nsites,natom*nsites),csizev(natom*nsites))
    !At the beginning each distance is a separate class and only
    !pairs with at least one atom in the primitive cell are considered
    allocate(myclass(0:nsites-1,0:nsites-1))
    nclass=0
    do ip=0,natom-1
       do is=ip,natom-1
          nclass=nclass+1
          if(is/=ip)then 
             csizev(nclass)=2
             patom(1,nclass)=ip; satom(1,nclass)=is
             patom(2,nclass)=is; satom(2,nclass)=ip
             myclass(ip,is)=nclass; myclass(is,ip)=nclass
          else 
             csizev(nclass)=1 
             patom(1,nclass)=ip; satom(1,nclass)=is
             myclass(ip,is)=nclass
          endif
       enddo
       do is=natom,nsites-1
          nclass=nclass+1
          patom(1,nclass)=ip; satom(1,nclass)=is
          csizev(nclass)=1
          myclass(ip,is)=nclass
       enddo
    enddo
    !Loop over symmetry operations. The "+1" symm op is pair permutation.
    do isymm=1,nsymm+1
       !Loop over classes (for the first symm op, classes are made
       !of all individual atom pairs in which the first atom lies
       !in the primitive cell and the second anywhere inside the supercell)
       do iclass=1,nclass
          istart=1
          do 
             csize=csizev(iclass)
             csizenew=csize
             do id=istart,csize
                !Map the atoms in the class under the symm operation
                if(isymm==nsymm+1)then
                   ip=satom(id,iclass)
                   is=patom(id,iclass)
                else
                   ip=SymmOp%map_symm(patom(id,iclass),isymm)
                   is=SymmOp%map_symm(satom(id,iclass),isymm)
                endif
                !translback is a vector that contain which translation is necessary
                !to bring the site back to its position inside the primitive cell
                itransl=SymmOp%translback(ip)
                ip2=SymmOp%translate(ip,itransl)
                is2=SymmOp%translate(is,itransl)
                jclass=myclass(ip2,is2)
                if(jclass/=iclass)then
                   jstart=csizenew
                   csizenew=csizenew+csizev(jclass)
                   do idj=1,csizev(jclass)
                      ip=patom(idj,jclass); is=satom(idj,jclass)
                      myclass(ip,is)=iclass
                      patom(jstart+idj,iclass)=ip
                      satom(jstart+idj,iclass)=is
                   enddo
                   csizev(jclass)=0
                endif
             enddo
             if(csizenew==csize)exit
             istart=csizev(iclass)+1
             csizev(iclass)=csizenew
          enddo
       enddo
    enddo
    mclass=0
    do i=1,nclass
       if(csizev(i)>0)mclass=mclass+1
       do j=1,csizev(i)
          ip=patom(j,i)
          is=satom(j,i)
          myclass(ip,is)=mclass
          do it=1,ntransl-1
             ip_transl=SymmOp%translate(ip,it)
             is_transl=SymmOp%translate(is,it)
             myclass(ip_transl,is_transl)=mclass
          enddo
       enddo
    enddo
    deallocate(patom,satom,csizev)
    nclass=mclass
    allocate(lattice%class_label(nclass,4),lattice%class_size(nclass))
    do iclass=1,nclass
       prim:do i=0,natom-1
          super:do j=0,nsites-1
             if(myclass(i,j)==iclass)then
                lattice%class_label(iclass,1:3)=lattice%cartpos(1:3,j)-lattice%cartpos(1:3,i)
                lattice%class_label(iclass,4)=dble(i)
                exit prim
             endif
          enddo super
       enddo prim
    enddo
    lattice%class_size(:)=0
    do i=0,nsites-1
       do j=0,nsites-1
          iclass=myclass(i,j)
          lattice%class_size(iclass)=lattice%class_size(iclass)+1
       enddo
    enddo
    lattice%nclass=nclass
    lattice%myclass=>myclass
    lattice%analyzed=.true.
  end subroutine construct_lattice_classes



  !----------------------------------------------------------------------------------------
  ! Create classes of equivalent k-points.
  !----------------------------------------------------------------------------------------
  subroutine construct_recip_lattice_classes(SymmOp,recip_lattice)
    use DQMC_RECLATT
    integer               :: i,ip,j,isymm,iclass,istart,csize,csizenew,&
         &jclass,idj,id,mclass,ip2,jstart,&
         &nsymm,nkpts
    integer,allocatable   :: patom(:,:),csizev(:)
    integer, pointer      :: myclass_k(:)
    type(symm_operations) :: SymmOp
    type(recip_lattice_t) :: recip_lattice

    if(.not.SymmOp%recip_lattice_mapped)stop'Need to map symmetries over lattice before classes (reciprocal)'

    !initialize local variables
    nsymm=SymmOp%nsymm
    nkpts=recip_lattice%nkpts
    if(SymmOp%TimeRev)nsymm=nsymm+1

    !Classes in k-space
    allocate(csizev(nkpts),myclass_k(nkpts),patom(nkpts,nkpts))
    csizev(:)=1
    do ip=1,nkpts
       csizev(ip)=1
       myclass_k(ip)=ip
       patom(1,ip)=ip
    enddo
    do isymm=1,nsymm
       do iclass=1,nkpts
          istart=1
          do 
             csize=csizev(iclass)
             csizenew=csize
             do id=istart,csize
                ip2=SymmOp%map_symm_k(patom(id,iclass),isymm)
                jclass=myclass_k(ip2)
                if(jclass/=iclass)then
                   jstart=csizenew
                   csizenew=csizenew+csizev(jclass)
                   do idj=1,csizev(jclass)
                      ip=patom(idj,jclass)
                      myclass_k(ip)=iclass
                      patom(jstart+idj,iclass)=ip
                   enddo
                   csizev(jclass)=0
                endif
             enddo
             if(csizenew==csize)exit
             istart=csizev(iclass)+1
             csizev(iclass)=csizenew
          enddo
       enddo
    enddo
    mclass=0
    do i=1,nkpts
       if(csizev(i)>0)mclass=mclass+1
       do j=1,csizev(i)
          ip=patom(j,i)
          myclass_k(ip)=mclass
       enddo
    enddo
    allocate(recip_lattice%class_size_k(mclass),recip_lattice%class_repr_k(mclass))
    j=0
    do i=1,nkpts
       if(csizev(i)>0)then
          j=j+1
          recip_lattice%class_repr_k(j)=i
          recip_lattice%class_size_k(j)=csizev(i)
       endif
    enddo
    deallocate(patom,csizev)
    recip_lattice%nclass_k=mclass
    recip_lattice%myclass_k=>myclass_k
    recip_lattice%analyzed=.true.
  end subroutine construct_recip_lattice_classes




  !---------------------------------------------------------------------------------
  ! This routines construct my_class_b(ib,jb) where ib and jb are two bonds.
  ! my_class contains the symmetry class of the pair (ib,jb)
  !---------------------------------------------------------------------------------
  subroutine construct_bond_classes(Bonds,SymmOp)
    use DQMC_BONDS
    type(symm_operations),intent(in) :: SymmOp
    type(bonds_t),intent(inout) :: Bonds
    integer :: ib,nclass,ntotbond,ntotbondsq,isymm,iclass,istart,csize,csizenew,&
         & id,bx,by,jclass,jstart,idj,mclass,jb,i,j
    integer,pointer :: myclass(:,:)
    integer, allocatable :: bond1(:,:),bond2(:,:),csizev(:)

    if(.not.SymmOp%bonds_mapped)stop'Need to map bonds before analyzing symmetry'

    ntotbond=size(SymmOp%map_symm_b,1)
    ntotbondsq=(ntotbond**2+ntotbond)/2
    allocate(myclass(ntotbond,ntotbond),bond1(ntotbondsq,ntotbondsq), &
         &        bond2(ntotbondsq,ntotbondsq),csizev(ntotbondsq))
    !Initially Define classes as if all bonds were different
    nclass=0
    do ib=1,ntotbond
       nclass=nclass+1
       myclass(ib,ib)=nclass
       bond1(1,nclass)=ib
       bond2(1,nclass)=ib
       csizev(nclass)=1
    enddo
    do ib=1,ntotbond
       do jb=ib+1,ntotbond
          nclass=nclass+1
          myclass(ib,jb)=nclass
          myclass(jb,ib)=nclass
          bond1(1,nclass)=ib
          bond2(1,nclass)=jb
          bond1(2,nclass)=jb
          bond2(2,nclass)=ib
          csizev(nclass)=2
       enddo
    enddo
    !Try all symmetry operations
    do isymm=1,SymmOp%nsymm
       !on all classes
       do iclass=1,nclass
          istart=1
          !we now loop over the elements of a class.
          !This number is increased as we found new equivalent
          !elements. That's why the loop is split in !1! and !2!
          do  !1!
             csize=csizev(iclass)
             csizenew=csize
             do id=istart,csize !2!
                bx=SymmOp%map_symm_b(bond1(id,iclass),isymm)
                by=SymmOp%map_symm_b(bond2(id,iclass),isymm)
                jclass=myclass(bx,by)
                if(jclass/=iclass)then
                   jstart=csizenew
                   csizenew=csizenew+csizev(jclass)
                   do idj=1,csizev(jclass)
                      bx=bond1(idj,jclass); by=bond2(idj,jclass)
                      myclass(bx,by)=iclass
                      bond1(jstart+idj,iclass)=bx
                      bond2(jstart+idj,iclass)=by
                   enddo
                   csizev(jclass)=0
                endif
             enddo
             if(csizenew==csize)exit
             istart=csizev(iclass)+1
             csizev(iclass)=csizenew
          enddo
       enddo
    enddo
    mclass=0
    do i=1,nclass
       if(csizev(i)>0)mclass=mclass+1
       do j=1,csizev(i)
          bx=bond1(j,i)
          by=bond2(j,i)
          myclass(bx,by)=mclass
       enddo
    enddo
    Bonds%nclass_b=mclass
    Bonds%myclass_b=>myclass
    allocate(Bonds%class_size_b(mclass))
    Bonds%class_size_b(:)=0.d0
    do ib=1,ntotbond
       do jb=1,ntotbond
          Bonds%class_size_b(myclass(ib,jb))=Bonds%class_size_b(myclass(ib,jb))+1
       enddo
    enddo
    !write(*,*)
    !write(*,*)'BOND PAIRS AND CLASS'
    !do ib=1,ntotbond
    ! do jb=1,ntotbond
    !  write(*,*) ib,jb,'-->',myclass(ib,jb)
    ! enddo
    !enddo
    deallocate(bond1,bond2,csizev)
    Bonds%analyzed=.true.
  end subroutine construct_bond_classes




  !---------------------------------------------------------------------------------
  ! Construct classes for pairs
  !---------------------------------------------------------------------------------
  subroutine construct_pair_classes(Bonds,Pairs)
    use DQMC_BONDS
    type(bonds_t),intent(in)    :: Bonds
    type(pairing),intent(inout) :: Pairs
    integer  :: nc, np, ip, jp, ib, jb, ic, jc
    integer, allocatable :: bclass(:)
    nc=0
    np=Pairs%nbond
    allocate(Pairs%myclass_p(np,np),bclass(Bonds%nclass_b))
    bclass=0
    do ip=1,np
       ib=Pairs%bond_map(ip)
       do jp=ip,np
          jb=Pairs%bond_map(jp)
          ic=Bonds%myclass_b(ib,jb)
          jc=bclass(ic)
          if(jc==0)then
             nc=nc+1
             Pairs%myclass_p(ip,jp)=nc
             Pairs%myclass_p(jp,ip)=nc
             bclass(ic)=nc
          else
             Pairs%myclass_p(ip,jp)=jc
             Pairs%myclass_p(jp,ip)=jc
          endif
       enddo
    enddo
    Pairs%nclass_p=nc
    deallocate(bclass)
    allocate(Pairs%class_size_p(nc))
    Pairs%class_size_p=0
    do ip=1,np
       do jp=1,np
          jc = Pairs%myclass_p(ip,jp)
          Pairs%class_size_p(jc)=Pairs%class_size_p(jc)+1
       enddo
    enddo
  end subroutine construct_pair_classes



end module DQMC_SYMM
